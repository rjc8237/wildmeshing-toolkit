#include "AdaptiveTessellation.h"
#include <fastenvelope/FastEnvelope.h>
#include <igl/Timer.h>
#include <igl/predicates/predicates.h>
#include <igl/read_triangle_mesh.h>
#include <igl/remove_unreferenced.h>
#include <igl/writeDMAT.h>
#include <igl/write_triangle_mesh.h>
#include <lagrange/utils/timing.h>
#include <tbb/concurrent_vector.h>
#include <wmtk/utils/AMIPS2D.h>
#include <wmtk/utils/AMIPS2D_autodiff.h>
#include <Eigen/Core>
#include <lean_vtk.hpp>
#include <map>
#include <tracy/Tracy.hpp>
#include <wmtk/utils/TriQualityUtils.hpp>
#include <wmtk/utils/TupleUtils.hpp>
#include "global_intersection.h"
using namespace wmtk;
namespace adaptive_tessellation {
// TODO change this to accomodate new error
auto avg_edge_len = [](auto& m) {
    double avg_len = 0.0;
    auto edges = m.get_edges();
    for (auto& e : edges) avg_len += std::sqrt(m.get_length3d(e));
    return avg_len / edges.size();
};

void AdaptiveTessellation::set_parameters(
    const double target_accuracy,
    const double target_edge_length,
    const wmtk::Image& image,
    const WrappingMode wrapping_mode,
    const SAMPLING_MODE sampling_mode,
    const DISPLACEMENT_MODE displacement_mode,
    const ENERGY_TYPE energy_type,
    const EDGE_LEN_TYPE edge_len_type,
    const bool boundary_parameter)
{
    if (mesh_parameters.m_edge_length_type == EDGE_LEN_TYPE::ACCURACY ||
        mesh_parameters.m_edge_length_type == EDGE_LEN_TYPE::AREA_ACCURACY)
        mesh_parameters.m_accuracy_threshold = target_accuracy;
    mesh_parameters.m_quality_threshold = target_edge_length;
    // setting needs to be in the order of image-> displacement-> energy-> edge_length
    // set the image first since it is used for displacement and energy setting
    set_image_function(image, wrapping_mode);
    mesh_parameters.m_sampling_mode = sampling_mode;
    set_displacement(displacement_mode);
    set_energy(energy_type);
    set_edge_length_measurement(edge_len_type);
    mesh_parameters.m_boundary_parameter = boundary_parameter;
} // namespace adaptive_tessellation

void AdaptiveTessellation::set_parameters(
    const double target_edge_length,
    const std::function<DScalar(const DScalar&, const DScalar&)>& displacement_function,
    const EDGE_LEN_TYPE edge_len_type,
    const ENERGY_TYPE energy_type,
    const bool boundary_parameter)
{
    mesh_parameters.m_target_l = target_edge_length;
    mesh_parameters.m_get_z = displacement_function;
    set_energy(
        energy_type); // set the displacement_function first since it is used for energy setting
    set_edge_length_measurement(edge_len_type); // set the default to get_legnth_3d
    mesh_parameters.m_boundary_parameter = boundary_parameter;
}
// using boundary parametrization,
// find the vertex that are the start and end of each cruve and set them as fixed
// use the vertex curve-id and assign edge curve-id
void AdaptiveTessellation::set_fixed()
{
    for (int curve_id = 0; curve_id < mesh_parameters.m_boundary.num_curves(); curve_id++) {
        if (mesh_parameters.m_boundary.is_periodic(curve_id)) {
            // periodic curve has no fixed vertex
            continue;
        }
        // set the first and last vertex as fixed
        auto uv_first = mesh_parameters.m_boundary.t_to_uv(curve_id, 0.);
        auto uv_last = mesh_parameters.m_boundary.t_to_uv(
            curve_id,
            mesh_parameters.m_boundary.upper_bound(curve_id));
        // find the closest points to the uv_first and uv_last
        double dist_first = std::numeric_limits<double>::infinity();
        double dist_last = std::numeric_limits<double>::infinity();
        size_t first_vid = -1;
        size_t last_vid = -1;
        for (auto& v : get_vertices()) {
            if ((vertex_attrs[v.vid(*this)].pos - uv_first).squaredNorm() < dist_first) {
                dist_first = (vertex_attrs[v.vid(*this)].pos - uv_first).squaredNorm();
                first_vid = v.vid(*this);
            }
            if ((vertex_attrs[v.vid(*this)].pos - uv_last).squaredNorm() < dist_last) {
                dist_last = (vertex_attrs[v.vid(*this)].pos - uv_last).squaredNorm();
                last_vid = v.vid(*this);
            }
        }
        assert(dist_first < 1e-8);
        assert(dist_last < 1e-8);
        vertex_attrs[first_vid].fixed = true;
        vertex_attrs[last_vid].fixed = true;
    }
}

/// @brief set the v as feature vertex, that is, it is a boundary vertex and the  incident boundary vertices are not colinear
/// @param v
void AdaptiveTessellation::set_feature(Tuple& v)
{
    // default fixed is false
    vertex_attrs[v.vid(*this)].fixed = false;
    // only set the feature if it is a boundary vertex
    assert(is_boundary_vertex(v));
    // get 3 vertices

    std::vector<size_t> incident_boundary_verts;
    auto one_ring = get_one_ring_edges_for_vertex(v);
    for (const auto& e : one_ring) {
        if (is_boundary_edge(e)) incident_boundary_verts.emplace_back(e.vid(*this));
    }
    // every boundary vertex has to have 2 neighbors that are boundary vertices
    assert(incident_boundary_verts.size() == 2);
    const auto& p1 = vertex_attrs[v.vid(*this)].pos;
    const auto& p2 = vertex_attrs[incident_boundary_verts[0]].pos;
    const auto& p3 = vertex_attrs[incident_boundary_verts[1]].pos;
    // check if they are co linear. set fixed = true if not
    double costheta = ((p2 - p1).stableNormalized()).dot((p3 - p1).stableNormalized()); // cos theta
    double theta = std::acos(std::clamp<double>(costheta, -1, 1));
    if (theta <= M_PI / 2) vertex_attrs[v.vid(*this)].fixed = true;
}

// assuming the m_triwild_displacement in mesh_parameter has been set
void AdaptiveTessellation::set_energy(const ENERGY_TYPE energy_type)
{
    std::unique_ptr<Energy> energy_ptr;
    switch (energy_type) {
    case ENERGY_TYPE::AMIPS: energy_ptr = std::make_unique<wmtk::AMIPS>(); break;
    case ENERGY_TYPE::SYMDI: energy_ptr = std::make_unique<wmtk::SymDi>(); break;
    case ENERGY_TYPE::EDGE_LENGTH:
        energy_ptr = std::make_unique<wmtk::EdgeLengthEnergy>(mesh_parameters.m_get_z);
        break;
    case ENERGY_TYPE::EDGE_QUADRATURE:
        energy_ptr = std::make_unique<wmtk::AccuracyEnergy>(mesh_parameters.m_displacement);
        break;
    case ENERGY_TYPE::AREA_QUADRATURE:
        energy_ptr = std::make_unique<wmtk::AreaAccuracyEnergy>(mesh_parameters.m_displacement);
        break;
    }

    mesh_parameters.m_energy = std::move(energy_ptr);
}

void AdaptiveTessellation::set_edge_length_measurement(const EDGE_LEN_TYPE edge_len_type)
{
    mesh_parameters.m_edge_length_type = edge_len_type;
    switch (edge_len_type) {
    case EDGE_LEN_TYPE::LINEAR2D:
        mesh_parameters.m_get_length = [&](const Tuple& edge_tuple) -> double {
            return this->get_length2d(edge_tuple);
        };
        break;
    case EDGE_LEN_TYPE::LINEAR3D:
        mesh_parameters.m_get_length = [&](const Tuple& edge_tuple) -> double {
            return this->get_length3d(edge_tuple);
        };
        break;
    case EDGE_LEN_TYPE::N_IMPLICIT_POINTS:
        mesh_parameters.m_get_length = [&](const Tuple& edge_tuple) -> double {
            return this->get_length_n_implicit_points(edge_tuple);
        };
        break;
    case EDGE_LEN_TYPE::PT_PER_PIXEL:
        mesh_parameters.m_get_length = [&](const Tuple& edge_tuple) -> double {
            return this->get_length_1ptperpixel(edge_tuple);
        };
        break;
    case EDGE_LEN_TYPE::MIPMAP:
        mesh_parameters.m_get_length = [&](const Tuple& edge_tuple) -> double {
            return this->get_length_mipmap(edge_tuple);
        };
        break;
    case EDGE_LEN_TYPE::ACCURACY:
        mesh_parameters.m_get_length = [&](const Tuple& edge_tuple) -> double {
            return this->get_edge_accuracy_error(edge_tuple);
        };
        break;
    case EDGE_LEN_TYPE::AREA_ACCURACY:
        mesh_parameters.m_get_length = [&](const Tuple& edge_tuple) -> double {
            return this->get_area_accuracy_error(edge_tuple);
        };
        break;
    }
}

void AdaptiveTessellation::set_image_function(
    const wmtk::Image& image,
    const WrappingMode wrapping_mode)
{
    mesh_parameters.m_wrapping_mode = wrapping_mode;
    mesh_parameters.m_image = image;
    mesh_parameters.m_get_z = [this](const DScalar& u, const DScalar& v) -> DScalar {
        return this->mesh_parameters.m_image.get(u, v);
    };
    mesh_parameters.m_image_get_coordinate =
        [this](const double& x, const double& y) -> std::pair<int, int> {
        auto [xx, yy] = this->mesh_parameters.m_image.get_pixel_index(x, y);
        return {
            this->mesh_parameters.m_image.get_coordinate(xx, this->mesh_parameters.m_wrapping_mode),
            this->mesh_parameters.m_image.get_coordinate(
                yy,
                this->mesh_parameters.m_wrapping_mode)};
    };
    mesh_parameters.m_mipmap = wmtk::MipMap(image);
    mesh_parameters.m_mipmap.set_wrapping_mode(wrapping_mode);
}

void AdaptiveTessellation::set_displacement(const DISPLACEMENT_MODE displacement_mode)
{
    // needs to be called after m_image is initiated
    // can be also set depending on a user parameter that initialize different Displacement type
    std::shared_ptr<Displacement> displacement_ptr;
    switch (displacement_mode) {
    case DISPLACEMENT_MODE::MESH_3D: {
        std::array<wmtk::Image, 6> position_normal_images;
        for (size_t i = 0; i < 2; i++) {
            std::filesystem::path path = mesh_parameters.m_position_normal_paths[i];
            wmtk::logger().info("======= path {} {}", i, path);
            wmtk::split_and_save_3channels(path);
            std::string directory = path.parent_path().string();
            std::string file = path.filename().string();
            std::filesystem::path path_r(directory + "/" + file + "_r.exr");
            std::filesystem::path path_g(directory + "/" + file + "_g.exr");
            std::filesystem::path path_b(directory + "/" + file + "_b.exr");
            position_normal_images[i * 3].load(
                path_r,
                mesh_parameters.m_wrapping_mode,
                mesh_parameters.m_wrapping_mode);
            position_normal_images[i * 3 + 1].load(
                path_g,
                mesh_parameters.m_wrapping_mode,
                mesh_parameters.m_wrapping_mode);
            position_normal_images[i * 3 + 2].load(
                path_b,
                mesh_parameters.m_wrapping_mode,
                mesh_parameters.m_wrapping_mode);
        }
        displacement_ptr = std::make_shared<DisplacementMesh>(
            mesh_parameters.m_image,
            position_normal_images,
            mesh_parameters.m_sampling_mode,
            mesh_parameters.m_scale,
            mesh_parameters.m_offset);
        break;
    }
    case DISPLACEMENT_MODE::PLANE:
        displacement_ptr = std::make_shared<DisplacementPlane>(
            mesh_parameters.m_image,
            mesh_parameters.m_sampling_mode);
        break;
    default: break;
    }
    mesh_parameters.m_displacement = displacement_ptr;
}

void AdaptiveTessellation::set_projection()
{
    struct Data
    {
        RowMatrix2<Index> E;
        RowMatrix2<Scalar> V;
        lagrange::bvh::EdgeAABBTree<RowMatrix2<Scalar>, RowMatrix2<Index>, 2> aabb;
    };
    auto data = std::make_shared<Data>();
    data->E = get_bnd_edge_matrix();
    RowMatrix2<Scalar> V_aabb = Eigen::MatrixXd::Zero(vert_capacity(), 2);
    for (int i = 0; i < vert_capacity(); ++i) {
        V_aabb.row(i) << vertex_attrs[i].pos[0], vertex_attrs[i].pos[1];
    }
    data->V = V_aabb;
    data->aabb =
        lagrange::bvh::EdgeAABBTree<RowMatrix2<Scalar>, RowMatrix2<Index>, 2>(data->V, data->E);
    std::function<Eigen::RowVector2d(const Eigen::RowVector2d& p)> projection =
        [data = std::move(data)](const Eigen::RowVector2d& p) -> Eigen::RowVector2d {
        uint64_t ind = 0;
        double distance = 0.0;
        Eigen::RowVector2d p_ret;
        assert(data != nullptr);
        (data->aabb).get_closest_point(p, ind, p_ret, distance);
        return p_ret;
    };

    mesh_parameters.m_get_closest_point = std::move(projection);
}

void AdaptiveTessellation::set_vertex_world_positions()
{
    for (const auto& v : get_vertices()) {
        const size_t v_id = v.vid(*this);
        vertex_attrs[v_id].pos_world = mesh_parameters.m_displacement->get(
            vertex_attrs[v_id].pos[0],
            vertex_attrs[v_id].pos[1]);
    }

    std::map<size_t, size_t> paired_vertices; // mapping from removed to remaining vertex

    // find seam vertices
    for (const auto& t : get_faces()) {
        const auto f_id = t.fid(*this);
        for (size_t e = 0; e < 3; ++e) {
            if (face_attrs[f_id].mirror_edges[e].has_value()) {
                Tuple e_tuple1 = tuple_from_edge(f_id, e);
                if (!e_tuple1.is_ccw(*this)) {
                    e_tuple1 = e_tuple1.switch_vertex(*this);
                }
                Tuple e_tuple2 = face_attrs[f_id].mirror_edges[e].value();
                if (!e_tuple2.is_ccw(*this)) {
                    e_tuple2 = e_tuple2.switch_vertex(*this);
                }
                const size_t v0 = e_tuple1.vid(*this);
                const size_t v1 = e_tuple1.switch_vertex(*this).vid(*this);
                const size_t v2 = e_tuple2.vid(*this);
                const size_t v3 = e_tuple2.switch_vertex(*this).vid(*this);
                if (v0 < v3) {
                    paired_vertices[v3] = v0;
                } else if (v3 < v0) {
                    paired_vertices[v0] = v3;
                }

                if (v1 < v2) {
                    paired_vertices[v2] = v1;
                } else if (v2 < v1) {
                    paired_vertices[v1] = v2;
                }
            }
        }
    }

    // make sure that all vertices are paired with the partner that has the lowest index
    for (const auto& [v0, _] : paired_vertices) {
        while (paired_vertices.count(paired_vertices[v0]) != 0) {
            paired_vertices[v0] = paired_vertices[paired_vertices[v0]];
        }
    }

    // collect all positions that belong to the same seam vertex
    std::map<size_t, std::vector<Eigen::Vector3d>> map_id_to_pos_vec;
    std::map<size_t, std::vector<size_t>> map_id_to_ids;
    for (const auto& [v0, v1] : paired_vertices) {
        if (map_id_to_pos_vec.count(v1) == 0) {
            map_id_to_pos_vec[v1] = {vertex_attrs[v1].pos_world};
            map_id_to_ids[v1].push_back(v1);
        }
        map_id_to_pos_vec[v1].push_back(vertex_attrs[v0].pos_world);
        map_id_to_ids[v1].push_back(v0);
    }

    // compute averate positions
    for (const auto& [v, pos_vec] : map_id_to_pos_vec) {
        Eigen::Vector3d p(0, 0, 0);
        for (const auto& pp : pos_vec) {
            p += pp;
        }
        p /= pos_vec.size();
        for (const auto& vv : map_id_to_ids[v]) {
            vertex_attrs[vv].pos_world = p;
        }
    }
}

void AdaptiveTessellation::create_paired_seam_mesh_with_offset(
    const std::filesystem::path input_mesh_path,
    Eigen::MatrixXd& UV,
    Eigen::MatrixXi& F)
{
    // load uv coordinates and connectivities
    auto mesh = lagrange::io::load_mesh<lagrange::SurfaceMesh32d>(input_mesh_path);
    triangulate_polygonal_facets(mesh);
    auto& uv_attr = mesh.get_indexed_attribute<double>(lagrange::AttributeName::texcoord);
    UV = lagrange::matrix_view(uv_attr.values());
    F = lagrange::reshaped_view(uv_attr.indices(), 3).cast<int>();
    assert(mesh.is_triangle_mesh());
    assert(mesh.get_num_facets() == F.rows());

    // load 3d coordinates and connectivities for computing the offset and scaling
    Eigen::MatrixXd V3d;
    Eigen::MatrixXi F3d;
    igl::read_triangle_mesh(input_mesh_path.string(), V3d, F3d);

    const Eigen::MatrixXd box_min = V3d.colwise().minCoeff();
    const Eigen::MatrixXd box_max = V3d.colwise().maxCoeff();
    double max_comp = (box_max - box_min).maxCoeff();
    Eigen::MatrixXd scene_offset = -box_min;
    Eigen::MatrixXd scene_extent = box_max - box_min;
    scene_offset.array() -= (scene_extent.array() - max_comp) * 0.5;
    mesh_parameters.m_scale = max_comp;
    mesh_parameters.m_offset = scene_offset;

    wmtk::logger().info("///// #v : {} {}", UV.rows(), UV.cols());
    wmtk::logger().info("///// #f : {} {}", F.rows(), F.cols());

    wmtk::TriMesh m_3d;
    std::vector<std::array<size_t, 3>> tris;
    for (auto f = 0; f < F3d.rows(); f++) {
        std::array<size_t, 3> tri = {(size_t)F3d(f, 0), (size_t)F3d(f, 1), (size_t)F3d(f, 2)};
        tris.emplace_back(tri);
    }
    m_3d.create_mesh(V3d.rows(), tris);
    create_mesh(UV, F);
    // loop through faces
    // for each edge, get fid on both side.
    // if it's boundary edge, do nothing
    // if not, index to the matrix for F_3d, and check the 3 edges
    // populate E0 and E1 for boundary construction
    std::map<std::pair<size_t, size_t>, std::pair<size_t, size_t>> seam_edges;
    Eigen::MatrixXi E0(F.rows() * 3, 2), E1(F.rows() * 3, 2);
    int seam_edge_cnt = 0;
    // for each vertex that's seam vertex, assign same color for mirror vertices
    // build the mapping from uv_index to color
    // and mapping from color to uv_index
    color_to_uv_indices.reserve(UV.rows());
    assert(uv_index_to_color.empty());
    int current_color = 0;
    for (auto fi = 0; fi < m_3d.tri_capacity(); ++fi) {
        for (auto lvi1 = 0; lvi1 < 3; ++lvi1) {
            auto lvi2 = (lvi1 + 1) % 3;
            auto local_eid = 3 - lvi1 - lvi2;
            // construct the edge tuple of current vertex in 3d mesh
            Tuple edge1_3d = Tuple(F3d(fi, lvi1), local_eid, fi, m_3d);
            assert(F3d(fi, lvi1) == edge1_3d.vid(m_3d));
            if (!edge1_3d.switch_face(m_3d).has_value()) {
                // Boundary edge, skipping...
                continue;
            } else {
                auto edge2_3d = edge1_3d.switch_face(m_3d).value();
                auto fj = edge2_3d.fid(m_3d);
                size_t lvj1, lvj2;
                for (auto i = 0; i < 3; i++) {
                    if (F3d(fj, i) == edge1_3d.vid(m_3d)) lvj1 = i;
                    if (F3d(fj, i) == edge1_3d.switch_vertex(m_3d).vid(m_3d)) lvj2 = i;
                }

                assert(F3d(fi, lvi1) == F3d(fj, lvj1));
                assert(F3d(fi, lvi2) == F3d(fj, lvj2));
                // set up seam edge
                if ((F(fi, lvi1) != F(fj, lvj1)) || (F(fi, lvi2) != F(fj, lvj2))) {
                    // this is a seam. init the mirror_edge tuple
                    // the orientation of the mirror edges is inccw (half edge conventions)
                    // However, the edge tuple in operations have arbitraty orientation
                    // !!! need to check orientations of mirror edge in operations !!!!!
                    TriMesh::Tuple seam_edge_fj(F(fj, lvj2), (3 - lvj1 - lvj2), fj, *this);
                    if (!seam_edge_fj.is_ccw(*this))
                        seam_edge_fj = seam_edge_fj.switch_vertex(*this);
                    face_attrs[fi].mirror_edges[local_eid] =
                        std::make_optional<wmtk::TriMesh::Tuple>(seam_edge_fj);

                    TriMesh::Tuple seam_edge_fi(F(fi, lvi1), local_eid, fi, *this);
                    if (!seam_edge_fi.is_ccw(*this))
                        seam_edge_fi = seam_edge_fi.switch_vertex(*this);
                    face_attrs[fj].mirror_edges[(3 - lvj1 - lvj2)] =
                        std::make_optional<wmtk::TriMesh::Tuple>(seam_edge_fi);
                    // check if seam_edge_fj is already in the map
                    size_t ei0 = std::min(F(fi, lvi1), F(fi, lvi2));
                    size_t ei1 = std::max(F(fi, lvi1), F(fi, lvi2));
                    size_t ej0 = std::min(F(fj, lvj1), F(fj, lvj2));
                    size_t ej1 = std::max(F(fj, lvj1), F(fj, lvj2));
                    if (seam_edges.find(std::pair{ei0, ei1}) != seam_edges.end()) {
                        assert((seam_edges[{ei0, ei1}] == std::pair<size_t, size_t>(ej0, ej1)));
                    } else {
                        // if not add it to the map and add to the edge matrix
                        seam_edges[{ej0, ej1}] = std::pair<size_t, size_t>(ei0, ei1);
                        E0.row(seam_edge_cnt) << ei0, ei1;
                        E1.row(seam_edge_cnt) << ej0, ej1;
                        seam_edge_cnt++;
                    }
                }
            }
            // set up fixed t-junction vertices
            // set up t-junction vertex coloring
            // get the one ring edges of this tuple
            if (uv_index_to_color.find(F(fi, lvi1)) != uv_index_to_color.end()) {
                // already colored, skipping...
                continue;
            } else {
                auto edge2_3d = edge1_3d.switch_face(m_3d).value();
                auto fj_3d = edge2_3d.fid(m_3d);
                size_t lvj1_3d, lvj2_3d;
                for (auto i = 0; i < 3; i++) {
                    if (F3d(fj_3d, i) == edge1_3d.vid(m_3d)) lvj1_3d = i;
                    if (F3d(fj_3d, i) == edge1_3d.switch_vertex(m_3d).vid(m_3d)) lvj2_3d = i;
                }

                assert(F3d(fi, lvi1) == F3d(fj_3d, lvj1_3d));
                assert(F3d(fi, lvi2) == F3d(fj_3d, lvj2_3d));
                // edge1_3d is a seam edge
                if ((F(fi, lvi1) != F(fj_3d, lvj1_3d)) || (F(fi, lvi2) != F(fj_3d, lvj2_3d))) {
                    uv_index_to_color.insert({F(fi, lvi1), current_color});
                    if (current_color < color_to_uv_indices.size() &&
                        std::find(
                            color_to_uv_indices[current_color].begin(),
                            color_to_uv_indices[current_color].end(),
                            F(fi, lvi1)) == color_to_uv_indices[current_color].end())
                        color_to_uv_indices[current_color].emplace_back(F(fi, lvi1));
                    else
                        color_to_uv_indices.emplace_back(1, F(fi, lvi1));
                    for (auto& e_3d : m_3d.get_one_ring_edges_for_vertex(edge1_3d)) {
                        if (!e_3d.switch_face(m_3d).has_value()) {
                            // Boundary edge, skipping...
                            continue;
                        } else {
                            auto e2_3d = e_3d.switch_face(m_3d).value();
                            auto fj = e2_3d.fid(m_3d);
                            size_t lvj1, lvj2;
                            for (auto i = 0; i < 3; i++) {
                                if (F3d(fj, i) == edge1_3d.vid(m_3d)) lvj1 = i;
                                if (F3d(fj, i) == e2_3d.vid(m_3d)) lvj2 = i;
                            }
                            assert(F3d(fi, lvi1) == F3d(fj, lvj1));
                            // this is a mirror vertex at a seam edge
                            if (F(fi, lvi1) != F(fj, lvj1)) {
                                // the mirror vertex has not been colored yet
                                if (uv_index_to_color.find(F(fj, lvj1)) ==
                                    uv_index_to_color.end()) {
                                    // add the color or the primary vertex to the mirror vertex
                                    uv_index_to_color.insert(
                                        {F(fj, lvj1), uv_index_to_color[F(fi, lvi1)]});
                                }
                                // if the mirror vertex is not in the color busket, add it
                                if (std::find(
                                        color_to_uv_indices[uv_index_to_color[F(fi, lvi1)]].begin(),
                                        color_to_uv_indices[uv_index_to_color[F(fi, lvi1)]].end(),
                                        F(fj, lvj1)) ==
                                    color_to_uv_indices[uv_index_to_color[F(fi, lvi1)]].end())
                                    color_to_uv_indices[uv_index_to_color[F(fi, lvi1)]]
                                        .emplace_back(F(fj, lvj1));
                                assert(
                                    uv_index_to_color[F(fj, lvj1)] ==
                                    uv_index_to_color[F(fi, lvi1)]);
                            }
                        }
                    }
                    current_color++;
                }
            }
        }
    }
    color_to_uv_indices.resize(current_color);
    E0.conservativeResize(seam_edge_cnt, 2);
    E1.conservativeResize(seam_edge_cnt, 2);
    assert(E0.rows() == E1.rows());
    mesh_construct_boundaries(UV, F, E0, E1);
}

bool AdaptiveTessellation::invariants(const std::vector<Tuple>& new_tris)
{
    TriMesh::invariants(new_tris);
    if (mesh_parameters.m_has_envelope) {
        for (auto& t : new_tris) {
            std::array<Eigen::Vector3d, 3> tris;
            auto vs = oriented_tri_vertices(t);
            for (auto j = 0; j < 3; j++) {
                tris[j] << vertex_attrs[vs[j].vid(*this)].pos(0),
                    vertex_attrs[vs[j].vid(*this)].pos(1), 0.0;
            }
            if (mesh_parameters.m_envelope.is_outside(tris)) {
                wmtk::logger().debug("envelop invariant fail");
                return false;
            }
        }
    }

    for (auto& t : new_tris) {
        Eigen::Vector2d a, b, c;
        auto verts = oriented_tri_vertices(t);
        assert(verts.size() == 3);
        a = vertex_attrs[verts[0].vid(*this)].pos;
        b = vertex_attrs[verts[1].vid(*this)].pos;
        c = vertex_attrs[verts[2].vid(*this)].pos;

        // check both inverted and exact colinear
        if (wmtk::orient2d_t(a, b, c) != 1) {
            wmtk::logger().debug(
                "----{} false in orientation and collinear {}",
                t.fid(*this),
                wmtk::orient2d_t(a, b, c));
            wmtk::logger().debug("{} {} {}", a, b, c);
            // write_ply("rejected_split_" + std::to_string(t.fid(*this)) + ".ply");
            return false;
        }
    }
    return true;
}
std::vector<AdaptiveTessellation::Tuple> AdaptiveTessellation::new_edges_after(
    const std::vector<AdaptiveTessellation::Tuple>& tris) const
{
    std::vector<AdaptiveTessellation::Tuple> new_edges;

    for (const auto& t : tris) {
        for (auto j = 0; j < 3; j++) {
            new_edges.push_back(tuple_from_edge(t.fid(*this), j));
        }
    }
    wmtk::unique_edge_tuples(*this, new_edges);
    return new_edges;
}

void AdaptiveTessellation::create_mesh(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F)
{
    std::vector<Eigen::Vector3d> V_env;
    V_env.resize(V.rows());
    std::vector<Eigen::Vector3i> F_env;
    F_env.resize(F.rows());
    // Register attributes
    p_vertex_attrs = &vertex_attrs;
    p_face_attrs = &face_attrs;
    p_edge_attrs = &edge_attrs;
    // Convert from eigen to internal representation (TODO: move to utils and remove it from all
    // app)
    std::vector<std::array<size_t, 3>> tri(F.rows());
    for (int i = 0; i < F.rows(); i++) {
        F_env[i] << (size_t)F(i, 0), (size_t)F(i, 1), (size_t)F(i, 2);
        tri[i][0] = (size_t)F(i, 0);
        tri[i][1] = (size_t)F(i, 1);
        tri[i][2] = (size_t)F(i, 2);
    }
    // Initialize the trimesh class which handles connectivity
    wmtk::TriMesh::create_mesh(V.rows(), tri);
    // Save the vertex position in the vertex attributes
    for (unsigned i = 0; i < V.rows(); ++i) {
        vertex_attrs[i].pos << V.row(i)[0], V.row(i)[1];
        V_env[i] << V.row(i)[0], V.row(i)[1], 0.0;
    }
    for (const auto& tri : this->get_faces()) {
        assert(!is_inverted(tri));
    }
}

// set fixed vertices due to boundary
// set the curve-id for each edge
void AdaptiveTessellation::mesh_construct_boundaries(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& F,
    const Eigen::MatrixXi& E0,
    const Eigen::MatrixXi& E1)
{
    // construct the boundary map for boundary parametrization
    if (mesh_parameters.m_boundary_parameter)
        mesh_parameters.m_boundary.construct_boundaries(V, F, E0, E1);
    // mark boundary vertices as boundary_vertex
    // but this is not indiscriminatively rejected for all operations
    // other operations are conditioned on whether m_bnd_freeze is turned on
    // also obtain the boudnary parametrizatin t for each vertex
    // for now keep the per vertex curve-id. but this is now a edge property
    for (auto v : this->get_vertices()) {
        if (is_boundary_vertex(v)) {
            vertex_attrs[v.vid(*this)].boundary_vertex = is_boundary_vertex(v);
            set_feature(v);
            // one vertex can have more than one curve-id.
            // current curve-id ofr vertex is arbitrarily picked among them
            std::tie(vertex_attrs[v.vid(*this)].curve_id, vertex_attrs[v.vid(*this)].t) =
                mesh_parameters.m_boundary.uv_to_t(vertex_attrs[v.vid(*this)].pos);
        }
    }
    // after the boundary is constructed, set the start and end of each curve to be fixed
    set_fixed();
    // assign curve-id to each edge using the curve-it assigned for each vertex
    assign_edge_curveid();
    // TODO move this to the boundary unit test
    for (const auto& v : get_vertices()) {
        assert(vertex_attrs[v.vid(*this)].t >= 0);
        if (is_boundary_vertex(v))
            assert(
                (vertex_attrs[v.vid(*this)].pos - mesh_parameters.m_boundary.t_to_uv(
                                                      vertex_attrs[v.vid(*this)].curve_id,
                                                      vertex_attrs[v.vid(*this)].t))
                    .squaredNorm() < 1e-8);
    }
}

Eigen::Matrix<uint64_t, Eigen::Dynamic, 2, Eigen::RowMajor>
AdaptiveTessellation::get_bnd_edge_matrix()
{
    int num_bnd_edge = 0;
    for (const auto& e : get_edges()) {
        if (is_boundary_edge(e)) num_bnd_edge++;
    }
    Eigen::Matrix<uint64_t, Eigen::Dynamic, 2, Eigen::RowMajor> E(num_bnd_edge, 2);
    int i = 0;
    for (const auto& e : get_edges()) {
        if (is_boundary_edge(e)) {
            E.row(i) << (uint64_t)e.vid(*this), (uint64_t)e.switch_vertex(*this).vid(*this);
            i++;
        }
    }
    return E;
}

void AdaptiveTessellation::export_mesh(Eigen::MatrixXd& V, Eigen::MatrixXi& F) const
{
    V = Eigen::MatrixXd::Zero(vert_capacity(), 2);
    for (auto& t : get_vertices()) {
        auto i = t.vid(*this);
        V.row(i) = vertex_attrs[i].pos;
    }

    F = Eigen::MatrixXi::Constant(tri_capacity(), 3, -1);
    for (auto& t : get_faces()) {
        auto i = t.fid(*this);
        auto vs = oriented_tri_vertices(t);
        for (int j = 0; j < 3; j++) {
            F(i, j) = vs[j].vid(*this);
        }
    }
}

void AdaptiveTessellation::export_mesh_without_invalid_faces(Eigen::MatrixXd& V, Eigen::MatrixXi& F)
    const
{
    V = Eigen::MatrixXd::Zero(vert_capacity(), 2);
    for (auto& t : get_vertices()) {
        auto i = t.vid(*this);
        V.row(i) = vertex_attrs[i].pos;
    }

    const auto faces = get_faces();
    F = Eigen::MatrixXi::Constant(faces.size(), 3, -1);
    for (size_t i = 0; i < faces.size(); ++i) {
        const auto vs = oriented_tri_vertices(faces[i]);
        for (size_t j = 0; j < 3; ++j) {
            F(i, j) = vs[j].vid(*this);
        }
    }
}

void AdaptiveTessellation::export_mesh_with_displacement(Eigen::MatrixXd& V, Eigen::MatrixXi& F)
    const
{
    export_mesh_without_invalid_faces(V, F);
    const size_t rows = V.rows();
    Eigen::MatrixXd vertices_displaced = Eigen::MatrixXd::Zero(rows, 3);
    for (size_t i = 0; i < rows; ++i) {
        vertices_displaced.row(i) = mesh_parameters.m_displacement->get(V(i, 0), V(i, 1));
    }

    V = vertices_displaced;
}

void AdaptiveTessellation::export_mesh(
    Eigen::MatrixXd& V,
    Eigen::MatrixXi& F,
    Eigen::MatrixXd& VT,
    Eigen::MatrixXi& FT) const
{
    export_mesh_without_invalid_faces(VT, FT);

    V.resize(VT.rows(), 3);
    for (int i = 0; i < VT.rows(); i++) {
        const double& u = VT(i, 0);
        const double& v = VT(i, 1);
        V.row(i) = mesh_parameters.m_displacement->get(u, v);
    }

    F = FT;
    remove_seams(V, F);

    // get rid of unreferenced vertices in both meshes
    Eigen::MatrixXd V_buf;
    Eigen::MatrixXi F_buf;
    Eigen::MatrixXi map_old_to_new_v_ids;
    igl::remove_unreferenced(V, F, V_buf, F_buf, map_old_to_new_v_ids);
    V = V_buf;
    F = F_buf;

    igl::remove_unreferenced(VT, FT, V_buf, F_buf, map_old_to_new_v_ids);
    VT = V_buf;
    FT = F_buf;
}

void AdaptiveTessellation::remove_seams(Eigen::MatrixXd& V, Eigen::MatrixXi& F) const
{
    std::map<size_t, size_t> paired_vertices; // mapping from removed to remaining vertex

    // find seam vertices
    for (const auto& t : get_faces()) {
        const auto f_id = t.fid(*this);
        for (size_t e = 0; e < 3; ++e) {
            if (face_attrs[f_id].mirror_edges[e].has_value()) {
                Tuple e_tuple1 = tuple_from_edge(f_id, e);
                if (!e_tuple1.is_ccw(*this)) {
                    e_tuple1 = e_tuple1.switch_vertex(*this);
                }
                Tuple e_tuple2 = face_attrs[f_id].mirror_edges[e].value();
                if (!e_tuple2.is_ccw(*this)) {
                    e_tuple2 = e_tuple2.switch_vertex(*this);
                }
                const size_t v0 = e_tuple1.vid(*this);
                const size_t v1 = e_tuple1.switch_vertex(*this).vid(*this);
                const size_t v2 = e_tuple2.vid(*this);
                const size_t v3 = e_tuple2.switch_vertex(*this).vid(*this);
                if (v0 < v3) {
                    paired_vertices[v3] = v0;
                } else if (v3 < v0) {
                    paired_vertices[v0] = v3;
                }

                if (v1 < v2) {
                    paired_vertices[v2] = v1;
                } else if (v2 < v1) {
                    paired_vertices[v1] = v2;
                }
            }
        }
    }

    // make sure that all vertices are paired with the partner that has the lowest index
    for (const auto& [v0, _] : paired_vertices) {
        while (paired_vertices.count(paired_vertices[v0]) != 0) {
            paired_vertices[v0] = paired_vertices[paired_vertices[v0]];
        }
    }

    // collect all positions that belong to the same seam vertex
    std::map<size_t, std::vector<Eigen::Vector3d>> map_id_to_pos_vec;
    for (const auto& [v0, v1] : paired_vertices) {
        if (map_id_to_pos_vec.count(v1) == 0) {
            map_id_to_pos_vec[v1] = {V.row(v1)};
        }
        map_id_to_pos_vec[v1].push_back(V.row(v0));
    }

    // compute average positions
    for (const auto& [v, pos_vec] : map_id_to_pos_vec) {
        Eigen::Vector3d p(0, 0, 0);
        for (const auto& pp : pos_vec) {
            p += pp;
        }
        p /= pos_vec.size();
        V.row(v) = p;
    }

    Eigen::MatrixXi NF;
    NF.resize(F.rows(), F.cols());
    for (size_t i = 0; i < NF.rows(); ++i) {
        for (size_t j = 0; j < NF.cols(); ++j) {
            size_t new_v_id = F(i, j);
            if (paired_vertices.count(new_v_id) != 0) {
                new_v_id = paired_vertices[new_v_id];
            }
            assert(paired_vertices.count(new_v_id) == 0);
            NF(i, j) = new_v_id;
        }
    }

    // overwrite F
    F = NF;
}

void AdaptiveTessellation::export_seamless_mesh_with_displacement(
    Eigen::MatrixXd& V,
    Eigen::MatrixXi& F) const
{
    export_mesh_with_displacement(V, F);
    remove_seams(V, F);
}

void AdaptiveTessellation::write_obj(const std::string& path)
{
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;

    export_mesh(V, F);

    Eigen::MatrixXd V3 = Eigen::MatrixXd::Zero(V.rows(), 3);
    V3.leftCols(2) = V;

    igl::writeOBJ(path, V3, F);
    wmtk::logger().info("============>> current edge length {}", avg_edge_len(*this));
}

void AdaptiveTessellation::write_ply(const std::string& path)
{
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;

    export_mesh(V, F);

    Eigen::MatrixXd V3 = Eigen::MatrixXd::Zero(V.rows(), 3);
    V3.leftCols(2) = V;

    igl::writePLY(path, V3, F);
}

void AdaptiveTessellation::write_vtk(const std::string& path)
{
    std::vector<double> points;
    std::vector<int> elements;
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;

    export_mesh(V, F);
    wmtk::logger().info("=== # vertices {}", V.rows());
    for (int i = 0; i < V.rows(); i++) {
        auto p = mesh_parameters.m_displacement->get(V(i, 0), V(i, 1));
        points.emplace_back(p(0));
        points.emplace_back(p(1));
        points.emplace_back(p(2));
    }
    assert(points.size() == 3 * V.rows());

    for (const auto& e : get_edges()) {
        elements.emplace_back(e.vid(*this));
        elements.emplace_back(e.switch_vertex(*this).vid(*this));
    }
    assert(elements.size() == 2 * get_edges().size());
    // vector<double> scalar_field = {0., 1., 2.};
    // vector<double> vector_field = points;
    const int dim = 3;
    const int cell_size = 2;
    leanvtk::VTUWriter writer;

    std::vector<double> scalar_field;
    for (const auto& e : get_edges()) {
        if (!e.is_valid(*this)) continue;
        auto cost = mesh_parameters.m_get_length(e);
        // Eigen::Matrix<double, 2, 1> pos1 = vertex_attrs[e.vid(*this)].pos;
        // Eigen::Matrix<double, 2, 1> pos2 =
        // vertex_attrs[e.switch_vertex(*this).vid(*this)].pos; Eigen::Matrix<double, 2, 1>
        // posnew = (pos1 + pos2) * 0.5; cost -=
        //     (mesh_parameters.m_displacement->get_error_per_edge(pos1, posnew) +
        //      mesh_parameters.m_displacement->get_error_per_edge(posnew, pos2));
        scalar_field.emplace_back(cost);
    }
    writer.add_cell_scalar_field("scalar_field", scalar_field);
    // writer.add_vector_field("vector_field", vector_field, dim);
    writer.write_surface_mesh(path, dim, cell_size, points, elements);
}
/// @brief write vtu with elements represent per face attributes
/// @param path
void AdaptiveTessellation::write_perface_vtk(const std::string& path)
{
    std::vector<double> points;
    std::vector<int> elements;
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;

    export_mesh(V, F);
    for (int i = 0; i < V.rows(); i++) {
        auto p = mesh_parameters.m_displacement->get(V(i, 0), V(i, 1));
        points.emplace_back(p(0));
        points.emplace_back(p(1));
        points.emplace_back(p(2));
    }
    assert(points.size() == 3 * V.rows());

    for (const auto& f : get_faces()) {
        auto its = oriented_tri_vids(f);
        elements.emplace_back(its[0]);
        elements.emplace_back(its[1]);
        elements.emplace_back(its[2]);
    }
    assert(elements.size() == 3 * get_faces().size());
    // vector<double> scalar_field = {0., 1., 2.};
    // vector<double> vector_field = points;
    const int dim = 3;
    const int cell_size = 3;
    leanvtk::VTUWriter writer;

    std::vector<double> scalar_field2;
    for (const auto& f : get_faces()) {
        if (!f.is_valid(*this)) continue;
        auto vids1 = oriented_tri_vertices(f);
        Eigen::Matrix<double, 3, 2, Eigen::RowMajor> triangle1;
        for (int i = 0; i < 3; i++) {
            triangle1.row(i) = vertex_attrs[vids1[i].vid(*this)].pos;
        }
        if (wmtk::polygon_signed_area(triangle1) < 0) {
            Eigen::Matrix<double, 1, 2, Eigen::RowMajor> tmp;
            tmp = triangle1.row(0);
            triangle1.row(0) = triangle1.row(1);
            triangle1.row(1) = tmp;
        }
        auto error = mesh_parameters.m_displacement->get_error_per_triangle(triangle1);

        scalar_field2.emplace_back(error);
    }
    writer.add_cell_scalar_field("scalar_field", scalar_field2);
    // writer.add_vector_field("vector_field", vector_field, dim);

    writer.write_surface_mesh(path, dim, cell_size, points, elements);
}

void AdaptiveTessellation::write_displaced_obj(
    const std::string& path,
    const std::function<double(double, double)>& displacement)
{
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;

    export_mesh(V, F);

    Eigen::MatrixXd V3 = Eigen::MatrixXd::Zero(V.rows(), 3);
    for (int i = 0; i < V.rows(); i++) {
        V3.row(i) << V(i, 0), V(i, 1), displacement(V(i, 0), V(i, 1));
    }

    igl::writeOBJ(path, V3, F);
    wmtk::logger().info("============>> current edge length {}", avg_edge_len(*this));
}

void AdaptiveTessellation::write_displaced_obj(
    const std::string& path,
    const std::shared_ptr<wmtk::Displacement> displacement)
{
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;

    export_mesh(V, F);
    auto rows = V.rows();
    Eigen::MatrixXd V3d = Eigen::MatrixXd::Zero(rows, 3);
    for (int i = 0; i < rows; i++) {
        V3d.row(i) = displacement->get(V(i, 0), V(i, 1));
        // wmtk::logger().info("progress: {}/{}", i, rows);
    }
    igl::writeOBJ(path, V3d, F);
    wmtk::logger().info("============>> current edge length {}", avg_edge_len(*this));
}

void AdaptiveTessellation::write_displaced_seamless_obj(
    const std::string& path,
    const std::shared_ptr<wmtk::Displacement> displacement)
{
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;

    export_seamless_mesh_with_displacement(V, F);

    igl::writeOBJ(path, V, F);
    wmtk::logger().info("============>> current edge length {}", avg_edge_len(*this));
}

void AdaptiveTessellation::write_world_obj(
    const std::string& path,
    const std::shared_ptr<wmtk::Displacement> displacement)
{
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;

    export_seamless_mesh_with_displacement(V, F);

    for (int i = 0; i < V.rows(); i++) {
        V.row(i) = vertex_attrs[i].pos_world;
    }

    igl::writeOBJ(path, V, F);
    wmtk::logger().info("============>> current edge length {}", avg_edge_len(*this));
}

void AdaptiveTessellation::write_obj_with_texture_coords(const std::string& path)
{
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    Eigen::MatrixXd CN;
    Eigen::MatrixXd FN;
    Eigen::MatrixXd VT;
    Eigen::MatrixXi FT;
    export_mesh(V, F, VT, FT);
    igl::writeOBJ(path, V, F, CN, FN, VT, FT);
    wmtk::logger().info("============>> current edge length {}", avg_edge_len(*this));
}

double AdaptiveTessellation::get_length2d(const Tuple& edge_tuple) const
{
    auto v1 = edge_tuple.vid(*this);
    auto v2 = edge_tuple.switch_vertex(*this).vid(*this);
    return (vertex_attrs[v1].pos - vertex_attrs[v2].pos).stableNorm();
}

double AdaptiveTessellation::get_length3d(const Tuple& edge_tuple) const
{
    auto v1 = edge_tuple.vid(*this);
    auto v2 = edge_tuple.switch_vertex(*this).vid(*this);
    auto v13d = mesh_parameters.m_project_to_3d(vertex_attrs[v1].pos(0), vertex_attrs[v1].pos(1));
    auto v23d = mesh_parameters.m_project_to_3d(vertex_attrs[v2].pos(0), vertex_attrs[v2].pos(1));
    return (v13d - v23d).stableNorm();
}

double AdaptiveTessellation::get_length_n_implicit_points(const Tuple& edge_tuple) const
{
    const auto& v1 = edge_tuple.vid(*this);
    const auto& v2 = edge_tuple.switch_vertex(*this).vid(*this);
    const auto& v12d = vertex_attrs[v1].pos;
    const auto& v22d = vertex_attrs[v2].pos;

    double length = 0.0;

    // add 3d displacement. add n implicit points to approximate quadrature of the curve
    std::vector<Eigen::Vector3d> quadrature;
    for (size_t i = 0; i < 6; i++) {
        auto tmp_v2d = v12d * (5 - i) / 5. + v22d * i / 5.;
        quadrature.emplace_back(mesh_parameters.m_project_to_3d(tmp_v2d(0), tmp_v2d(1)));
    }
    for (size_t i = 0; i < 5; i++) {
        length += (quadrature[i] - quadrature[i + 1]).stableNorm();
    }
    return length;
}

double AdaptiveTessellation::get_length_1ptperpixel(const Tuple& edge_tuple) const
{
    const auto& vid1 = edge_tuple.vid(*this);
    const auto& vid2 = edge_tuple.switch_vertex(*this).vid(*this);
    const auto& v12d = vertex_attrs[vid1].pos;
    const auto& v22d = vertex_attrs[vid2].pos;

    double length = 0.0;
    // get the pixel index of p1 and p2
    const auto [xx1, yy1] = mesh_parameters.m_image_get_coordinate(v12d(0), v12d(1));
    const auto [xx2, yy2] = mesh_parameters.m_image_get_coordinate(v22d(0), v22d(1));
    // get all the pixels in between p1 and p2
    const size_t pixel_num = std::max(abs(xx2 - xx1), abs(yy2 - yy1));
    // sum up the length
    // add 3d displacement. add n implicit points to approximate quadrature of the curve
    if (pixel_num <= 0) return -1.;
    assert(pixel_num > 0);
    std ::vector<Eigen::Vector3d> quadrature;
    for (size_t i = 0; i < (pixel_num + 1); i++) {
        const double u = static_cast<double>(i) / pixel_num;
        auto tmp_v2d = v12d * (1. - u) + v22d * u;
        quadrature.emplace_back(mesh_parameters.m_project_to_3d(tmp_v2d(0), tmp_v2d(1)));
    }
    for (size_t i = 0; i < pixel_num; i++) {
        length += (quadrature[i] - quadrature[i + 1]).stableNorm();
    }
    return length;
}

double AdaptiveTessellation::get_length_mipmap(const Tuple& edge_tuple) const
{
    const auto& vid1 = edge_tuple.vid(*this);
    const auto& vid2 = edge_tuple.switch_vertex(*this).vid(*this);
    const auto& v12d = vertex_attrs[vid1].pos;
    const auto& v22d = vertex_attrs[vid2].pos;
    int idx = 0;
    int pixel_num = -1;

    std::tie(idx, pixel_num) = mesh_parameters.m_mipmap.get_mipmap_level_pixelnum(v12d, v22d);
    const auto& image = mesh_parameters.m_mipmap.get_image(idx);
    double length = 0.0;

    if (pixel_num <= 0) return -1.;
    // sum up the length
    // add 3d displacement. add n implicit points to approximate quadrature of the curve
    std::vector<Eigen::Vector3d> quadrature;
    for (int i = 0; i < (pixel_num + 1); i++) {
        auto tmp_v2d = v12d * (pixel_num - i) / pixel_num + v22d * i / pixel_num;
        double z = image.get(tmp_v2d(0), tmp_v2d(1));
        quadrature.emplace_back(tmp_v2d(0), tmp_v2d(1), z);
    }
    for (size_t i = 0; i < pixel_num; i++) {
        length += (quadrature[i] - quadrature[i + 1]).stableNorm();
    }
    return length;
}

double AdaptiveTessellation::get_quality(const Tuple& loc, int idx) const
{
    // Global ids of the vertices of the triangle
    auto its = oriented_tri_vids(loc);

    // Temporary variable to store the stacked coordinates of the triangle
    std::array<double, 6> T;
    auto energy = -1.;
    for (size_t k = 0; k < 3; k++)
        for (size_t j = 0; j < 2; j++) T[k * 2 + j] = vertex_attrs[its[k]].pos[j];

    // Energy evaluation
    // energy = wmtk::AMIPS2D_energy(T);
    wmtk::State state = {};
    state.input_triangle = T;
    state.scaling = mesh_parameters.m_target_l;
    state.idx = idx;
    mesh_parameters.m_energy->eval(state);
    energy = state.value;

    // Filter for numerical issues
    if (std::isinf(energy) || std::isnan(energy)) return mesh_parameters.MAX_ENERGY;

    return energy;
}

double AdaptiveTessellation::get_edge_accuracy_error(const Tuple& edge_tuple) const
{
    const auto& vid1 = edge_tuple.vid(*this);
    const auto& vid2 = edge_tuple.switch_vertex(*this).vid(*this);
    const auto& v12d = vertex_attrs[vid1].pos;
    const auto& v22d = vertex_attrs[vid2].pos;

    return mesh_parameters.m_displacement->get_error_per_edge(v12d, v22d);
}

double AdaptiveTessellation::get_area_accuracy_error_per_face(const Tuple& edge_tuple) const
{
    auto vids = oriented_tri_vertices(edge_tuple);
    Eigen::Matrix<double, 3, 2, Eigen::RowMajor> triangle;
    for (int i = 0; i < 3; i++) {
        triangle.row(i) = vertex_attrs[vids[i].vid(*this)].pos;
    }
    auto triangle_area = wmtk::polygon_signed_area(triangle);
    if (triangle_area < 0.) {
        Eigen::Matrix<double, 1, 2, Eigen::RowMajor> tmp;
        tmp = triangle.row(0);
        triangle.row(0) = triangle.row(1);
        triangle.row(1) = tmp;
    }
    return mesh_parameters.m_displacement->get_error_per_triangle(triangle);
}

double AdaptiveTessellation::get_area_accuracy_error_per_face_triangle_matrix(
    Eigen::Matrix<double, 3, 2, Eigen::RowMajor> triangle) const
{
    auto triangle_area = wmtk::polygon_signed_area(triangle);
    if (triangle_area < 0.) {
        Eigen::Matrix<double, 1, 2, Eigen::RowMajor> tmp;
        tmp = triangle.row(0);
        triangle.row(0) = triangle.row(1);
        triangle.row(1) = tmp;
    }
    return mesh_parameters.m_displacement->get_error_per_triangle(triangle);
}

double AdaptiveTessellation::get_area_accuracy_error(const Tuple& edge_tuple) const
{
    double error = 0.;
    double error1 = get_area_accuracy_error_per_face(edge_tuple);
    double error2 = 0.;
    auto tmp_tuple = edge_tuple.switch_face(*this);
    if (tmp_tuple.has_value()) {
        error2 = get_area_accuracy_error_per_face(tmp_tuple.value());
    } else {
        if (is_seam_edge(edge_tuple))
            error2 = get_area_accuracy_error_per_face(get_oriented_mirror_edge(edge_tuple));
        else
            error2 = error1;
    }
    if (mesh_parameters.m_split_absolute_error_metric) {
        error = (error1 + error2) * get_length2d(edge_tuple);
    } else {
        double e_before = error1 + error2;
        Eigen::Matrix<double, 3, 2, Eigen::RowMajor> triangle;
        double e_after, error_after_1, error_after_2, error_after_3, error_after_4;

        triangle.row(0) = vertex_attrs[edge_tuple.vid(*this)].pos;
        triangle.row(1) = (vertex_attrs[edge_tuple.vid(*this)].pos +
                           vertex_attrs[edge_tuple.switch_vertex(*this).vid(*this)].pos) /
                          2.;
        triangle.row(2) =
            vertex_attrs[edge_tuple.switch_edge(*this).switch_vertex(*this).vid(*this)].pos;

        error_after_1 = get_area_accuracy_error_per_face_triangle_matrix(triangle);
        triangle.row(0) = vertex_attrs[edge_tuple.switch_vertex(*this).vid(*this)].pos;
        error_after_2 = get_area_accuracy_error_per_face_triangle_matrix(triangle);
        if (edge_tuple.switch_face(*this).has_value()) {
            triangle.row(2) = vertex_attrs[(edge_tuple.switch_face(*this).value())
                                               .switch_edge(*this)
                                               .switch_vertex(*this)
                                               .vid(*this)]
                                  .pos;
            error_after_3 = get_area_accuracy_error_per_face_triangle_matrix(triangle);
            triangle.row(0) = vertex_attrs[edge_tuple.vid(*this)].pos;
            error_after_4 = get_area_accuracy_error_per_face_triangle_matrix(triangle);
        } else {
            error_after_3 = error_after_1;
            error_after_4 = error_after_2;
        }
        e_after = error_after_1 + error_after_2 + error_after_3 + error_after_4;
        error = e_before - e_after;
    }
    return error;
}

void AdaptiveTessellation::get_nminfo_for_vertex(const Tuple& v, wmtk::NewtonMethodInfo& nminfo)
    const
{
    auto is_inverted_coordinates = [](auto& A, auto& B, auto& primary_pos) {
        auto res = igl::predicates::orient2d(A, B, primary_pos);
        if (res != igl::predicates::Orientation::POSITIVE)
            return true;
        else
            return false;
    };
    std::vector<Tuple> one_ring_tris = get_one_ring_tris_for_vertex(v);
    nminfo.curve_id = vertex_attrs[v.vid(*this)].curve_id;
    nminfo.target_length = mesh_parameters.m_target_l;
    nminfo.neighbors.resize(one_ring_tris.size(), 4);
    for (auto i = 0; i < one_ring_tris.size(); i++) {
        const Tuple& tri = one_ring_tris[i];
        assert(!is_inverted(tri));
        std::array<Tuple, 3> local_tuples = oriented_tri_vertices(tri);
        for (size_t j = 0; j < 3; j++) {
            if (local_tuples[j].vid(*this) == v.vid(*this)) {
                const Eigen::Vector2d& v2 = vertex_attrs[local_tuples[(j + 1) % 3].vid(*this)].pos;
                const Eigen::Vector2d& v3 = vertex_attrs[local_tuples[(j + 2) % 3].vid(*this)].pos;
                nminfo.neighbors.row(i) << v2(0), v2(1), v3(0), v3(1);
                assert(!is_inverted_coordinates(v2, v3, vertex_attrs[v.vid(*this)].pos));
                // sanity check. Should not be inverted
            }
        }
        assert(one_ring_tris.size() == nminfo.neighbors.rows());
    }
}

std::pair<double, Eigen::Vector2d> AdaptiveTessellation::get_one_ring_energy(const Tuple& v) const
{
    std::vector<wmtk::TriMesh::Tuple> one_ring_tris = get_one_ring_tris_for_vertex(v);
    assert(one_ring_tris.size() > 0);
    wmtk::DofVector dofx;
    if (is_boundary_vertex(v) && mesh_parameters.m_boundary_parameter) {
        dofx.resize(1);
        dofx[0] = vertex_attrs[v.vid(*this)].t; // t
    } else {
        dofx.resize(2); // uv;
        dofx = vertex_attrs[v.vid(*this)].pos;
    }

    // infomation needed for newton's method
    // multiple nminfo for seam vertices
    std::vector<wmtk::NewtonMethodInfo> nminfos;
    // push in current vertex's nminfo
    wmtk::NewtonMethodInfo primary_nminfo;
    get_nminfo_for_vertex(v, primary_nminfo);
    nminfos.emplace_back(primary_nminfo);
    // check if it is seam by getting the one ring edges
    // add the nminfo for each seam edge
    // keep a list of mirror vertices to update the dofx later
    std::vector<wmtk::TriMesh::Tuple> mirror_vertices;
    for (auto& e : get_one_ring_edges_for_vertex(v)) {
        if (is_seam_edge(e)) {
            wmtk::NewtonMethodInfo nminfo;
            assert(e.switch_vertex(*this).vid(*this) == v.vid(*this));
            wmtk::TriMesh::Tuple mirror_v = get_mirror_vertex(e.switch_vertex(*this));
            mirror_vertices.emplace_back(mirror_v);
            // collect the triangles for invariants check
            for (auto& mirror_v_tri : get_one_ring_tris_for_vertex(mirror_v)) {
                one_ring_tris.emplace_back(mirror_v_tri);
            }
            get_nminfo_for_vertex(mirror_v, nminfo);
            nminfos.push_back(nminfo);
        }
    }
    State state = {};
    state.dofx = dofx;
    wmtk::optimization_state_update(
        *mesh_parameters.m_energy,
        nminfos,
        mesh_parameters.m_boundary,
        state);

    return {state.value, state.gradient};
}

Eigen::VectorXd AdaptiveTessellation::get_quality_all_triangles()
{
    // Use a concurrent vector as for_each_face is parallel
    tbb::concurrent_vector<double> quality;
    quality.reserve(vertex_attrs.size());

    // Evaluate quality in parallel
    for_each_face([&](auto& f) { quality.push_back(get_quality(f)); });

    // Copy back in a VectorXd
    Eigen::VectorXd ret(quality.size());
    for (unsigned i = 0; i < quality.size(); ++i) ret[i] = quality[i];
    return ret;
}

bool AdaptiveTessellation::is_inverted(const Tuple& loc) const
{
    // Get the vertices ids
    auto vs = oriented_tri_vertices(loc);

    igl::predicates::exactinit();

    // Use igl for checking orientation
    auto res = igl::predicates::orient2d(
        vertex_attrs[vs[0].vid(*this)].pos,
        vertex_attrs[vs[1].vid(*this)].pos,
        vertex_attrs[vs[2].vid(*this)].pos);
    // The element is inverted if it not positive (i.e. it is negative or it is degenerate)
    return (res != igl::predicates::Orientation::POSITIVE);
}

void AdaptiveTessellation::mesh_improvement(int max_its)
{
    auto start_time = lagrange::get_timestamp();

    [[maybe_unused]] igl::Timer timer;
    [[maybe_unused]] double avg_len = 0.0;
    [[maybe_unused]] double pre_avg_len = 0.0;
    [[maybe_unused]] double pre_max_energy = -1.0;
    [[maybe_unused]] double old_average = 1e-4;
    wmtk::logger().info("target len {}", mesh_parameters.m_target_l);
    wmtk::logger().info("current length {}", avg_edge_len(*this));
    // mesh_parameters.js_log["edge_length_avg_start"] = avg_edge_len(*this);
    for (int it = 0; it < max_its; it++) {
        wmtk::logger().info("\n========it {}========", it);

        ///energy check
        wmtk::logger().info(
            "current max energy {} stop energy {}",
            mesh_parameters.m_max_energy,
            mesh_parameters.m_stop_energy);
        wmtk::logger().info("current length {}", avg_edge_len(*this));

        split_all_edges();
        assert(invariants(get_faces()));
        auto split_finish_time = lagrange::get_timestamp();
        mesh_parameters.js_log["iteration_" + std::to_string(it)]["split time"] =
            lagrange::timestamp_diff_in_seconds(start_time, split_finish_time);
        // consolidate_mesh();
        write_displaced_obj(
            mesh_parameters.m_output_folder + "/after_split_" + std::to_string(it) + ".obj",
            mesh_parameters.m_displacement);
        write_obj(
            mesh_parameters.m_output_folder + "/after_split_" + std::to_string(it) + "2d.obj");
        displace_self_intersection_free(*this);
        write_world_obj(
            mesh_parameters.m_output_folder + "/after_split_" + std::to_string(it) +
                "3d_intersection_free.obj",
            mesh_parameters.m_displacement);

        swap_all_edges();
        assert(invariants(get_faces()));
        auto swap_finish_time = lagrange::get_timestamp();
        mesh_parameters.js_log["iteration_" + std::to_string(it)]["swap time"] =
            lagrange::timestamp_diff_in_seconds(split_finish_time, swap_finish_time);
        consolidate_mesh();
        write_displaced_obj(
            mesh_parameters.m_output_folder + "/after_swap_" + std::to_string(it) + ".obj",
            mesh_parameters.m_displacement);
        write_obj(mesh_parameters.m_output_folder + "/after_swap_" + std::to_string(it) + "2d.obj");

        collapse_all_edges();
        assert(invariants(get_faces()));
        consolidate_mesh();
        auto collapse_finish_time = lagrange::get_timestamp();
        mesh_parameters.js_log["iteration_" + std::to_string(it)]["collapse time"] =
            lagrange::timestamp_diff_in_seconds(swap_finish_time, collapse_finish_time);

        write_displaced_obj(
            mesh_parameters.m_output_folder + "/after_collapse_" + std::to_string(it) + ".obj",
            mesh_parameters.m_displacement);
        write_obj(
            mesh_parameters.m_output_folder + "/after_collapse_" + std::to_string(it) + "2d.obj");

        smooth_all_vertices();
        assert(invariants(get_faces()));
        auto smooth_finish_time = lagrange::get_timestamp();
        mesh_parameters.js_log["iteration_" + std::to_string(it)]["smooth time"] =
            lagrange::timestamp_diff_in_seconds(collapse_finish_time, smooth_finish_time);
        consolidate_mesh();
        write_displaced_obj(
            mesh_parameters.m_output_folder + "/after_smooth_" + std::to_string(it) + ".obj",
            mesh_parameters.m_displacement);
        write_obj(
            mesh_parameters.m_output_folder + "/after_smooth_" + std::to_string(it) + "2d.obj");

        auto avg_grad = (mesh_parameters.m_gradient / vert_capacity()).stableNorm();

        wmtk::logger().info(
            "++++++++v {} t {} avg gradient {}++++++++",
            vert_capacity(),
            tri_capacity(),
            avg_grad);

        mesh_parameters.js_log["iteration_" + std::to_string(it)]["avg_grad"] = avg_grad;
        mesh_parameters.js_log["iteration_" + std::to_string(it)]["num_v"] = vert_capacity();
        mesh_parameters.js_log["iteration_" + std::to_string(it)]["num_f"] = tri_capacity();
        mesh_parameters.js_log["iteration_" + std::to_string(it)]["energy_max"] =
            mesh_parameters.m_max_energy;
        mesh_parameters.js_log["iteration_" + std::to_string(it)]["edge_len_avg"] =
            avg_edge_len(*this);
        mesh_parameters.js_log["iteration_" + std::to_string(it)]["edge_len_target"] =
            mesh_parameters.m_target_l;
        if (avg_grad < 1e-4) {
            wmtk::logger().info(
                "!!!avg grad is less than 1e-4 !!! energy doesn't improve anymore. early stop"
                "itr {}, avg length {} ",
                it,
                avg_len);
            break;
        }
        mesh_parameters.m_gradient = Eigen::Vector2d(0., 0.);
        avg_len = avg_edge_len(*this);
        mesh_parameters.js_log["iteration_" + std::to_string(it)]["edge_len_avg_final"] = avg_len;
        if (abs(avg_grad - old_average) < 1e-4) break;
        old_average = avg_grad;
        if (mesh_parameters.m_target_l <= 0 &&
            mesh_parameters.m_max_energy < mesh_parameters.m_stop_energy) {
            break;
        }
        pre_avg_len = avg_len;
        pre_max_energy = mesh_parameters.m_max_energy;
        consolidate_mesh();
        FrameMark;
    }

    wmtk::logger().info(
        "/////final: max energy {} , avg len {} ",
        mesh_parameters.m_max_energy,
        avg_len);
    consolidate_mesh();
}
void AdaptiveTessellation::flatten_dofs(Eigen::VectorXd& v_flat)
{
    auto verts = get_vertices();
    v_flat.resize(verts.size() * 2);
    for (const auto& v : verts) {
        if (is_boundary_vertex(v) && mesh_parameters.m_boundary_parameter) {
            v_flat(v.vid(*this) * 2) = vertex_attrs[v.vid(*this)].t;
            v_flat(v.vid(*this) * 2 + 1) = std::numeric_limits<double>::infinity();
        } else {
            v_flat(v.vid(*this) * 2) = vertex_attrs[v.vid(*this)].pos(0);
            v_flat(v.vid(*this) * 2 + 1) = vertex_attrs[v.vid(*this)].pos(1);
        }
    }
}

// get the energy defined by edge_length_energy over each face of the mesh
// assuming the vert_capacity() == get_vertices.size()
double AdaptiveTessellation::get_mesh_energy(const Eigen::VectorXd& v_flat)
{
    double total_energy = 0;
    Eigen::MatrixXd energy_matrix;
    energy_matrix.resize(get_faces().size(), 2);
    for (auto& face : get_faces()) {
        // wmtk::logger().info("getting energy on {} ", f_cnt++);
        auto verts = oriented_tri_vertices(face);
        Eigen::Matrix3d v_matrix;
        v_matrix.setZero(3, 3);
        for (int i = 0; i < 3; i++) {
            const auto& vert = verts[i];
            if (is_boundary_vertex(vert) && mesh_parameters.m_boundary_parameter) {
                auto uv = mesh_parameters.m_boundary.t_to_uv(
                    vertex_attrs[vert.vid(*this)].curve_id,
                    v_flat[vert.vid(*this) * 2]);
                v_matrix.row(i) = mesh_parameters.m_project_to_3d(uv(0), uv(1));
            } else {
                auto u = v_flat[vert.vid(*this) * 2];
                auto v = v_flat[vert.vid(*this) * 2 + 1];
                v_matrix.row(i) = mesh_parameters.m_project_to_3d(u, v);
            }
        }

        assert(mesh_parameters.m_target_l != 0);
        double tri_energy = 0.;
        auto BA = v_matrix.row(1) - v_matrix.row(0);
        auto CA = v_matrix.row(2) - v_matrix.row(0);
        auto BC = v_matrix.row(1) - v_matrix.row(2);
        tri_energy += pow(BA.squaredNorm() - pow(mesh_parameters.m_target_l, 2), 2);
        tri_energy += pow(BC.squaredNorm() - pow(mesh_parameters.m_target_l, 2), 2);
        tri_energy += pow(CA.squaredNorm() - pow(mesh_parameters.m_target_l, 2), 2);

        energy_matrix(face.fid(*this), 0) = tri_energy;
        double area = (BA.cross(CA)).squaredNorm();
        double A_hat = 0.5 * (std::sqrt(3) / 2) * 0.5 *
                       pow(mesh_parameters.m_target_l, 2); // this is arbitrary now
        assert(A_hat > 0);
        if (area <= 0) {
            tri_energy += std::numeric_limits<double>::infinity();
        }
        if (area < A_hat) {
            assert((area / A_hat) < 1.0);
            tri_energy += -(area - A_hat) * (area - A_hat) * log(area / A_hat);
        }
        total_energy += tri_energy;

        energy_matrix(face.fid(*this), 1) = tri_energy;
    }
    igl::writeDMAT("mesh_energy.dmat", energy_matrix);

    return total_energy;
}
/// debugging
void AdaptiveTessellation::gradient_debug(int max_its)
{
    split_all_edges();
    assert(invariants(get_faces()));
    consolidate_mesh();
    write_obj("before_smooth.obj");
    double old_average = 1e-4;
    for (int it = 0; it < max_its; it++) {
        wmtk::logger().info("\n========it {}========", it);
        wmtk::logger().info("current length {}", avg_edge_len(*this));

        smooth_all_vertices();
        assert(invariants(get_faces()));
        consolidate_mesh();
        write_displaced_obj(
            "smooth_" + std::to_string(it) + ".obj",
            mesh_parameters.m_displacement);

        wmtk::logger().info(
            "++++++++v {} t {} avg gradient {}++++++++",
            vert_capacity(),
            tri_capacity(),
            mesh_parameters.m_gradient / vert_capacity());
        auto avg_grad = (mesh_parameters.m_gradient / vert_capacity()).stableNorm();

        mesh_parameters.js_log["iteration_" + std::to_string(it)]["avg_grad"] = avg_grad;
        Eigen::VectorXd v_flat;
        flatten_dofs(v_flat);
        mesh_parameters.js_log["iteration_" + std::to_string(it)]["energy"] =
            get_mesh_energy(v_flat);
        if (abs(avg_grad - old_average) < 1e-5) break;
        old_average = avg_grad;
        mesh_parameters.m_gradient = Eigen::Vector2d(0., 0.);
    }
}

bool AdaptiveTessellation::is_seam_edge(const TriMesh::Tuple& t) const
{
    return face_attrs[t.fid(*this)].mirror_edges[t.local_eid(*this)].has_value();
}
bool AdaptiveTessellation::is_seam_vertex(const TriMesh::Tuple& t) const
{
    // it is seam vertex if it's incident to a seam edge
    for (auto& e : get_one_ring_edges_for_vertex(t)) {
        if (is_seam_edge(e)) return true;
    }
    return false;
}

void AdaptiveTessellation::set_mirror_edge_data(
    const TriMesh::Tuple& primary_t,
    const TriMesh::Tuple& mirror_edge)
{
    face_attrs[primary_t.fid(*this)].mirror_edges[primary_t.local_eid(*this)] =
        mirror_edge.is_ccw(*this) ? mirror_edge : mirror_edge.switch_vertex(*this);
}
std::optional<TriMesh::Tuple> AdaptiveTessellation::get_sibling_edge(const TriMesh::Tuple& t) const
{
    if (is_boundary_edge(t)) {
        if (is_seam_edge(t)) {
            return std::make_optional<TriMesh::Tuple>(get_oriented_mirror_edge(t));
        } else {
            return std::nullopt;
        }
    } else {
        assert(t.switch_face(*this).has_value());
        return std::make_optional<TriMesh::Tuple>(
            t.switch_face(*this).value().switch_vertex(*this));
        // return the sibling edge that's of opposite diretion
    }
}

// given a seam edge retrieve its mirror edge in opposite direction (half egde conventions )
TriMesh::Tuple AdaptiveTessellation::get_oriented_mirror_edge(const TriMesh::Tuple& t) const
{
    assert(is_seam_edge(t));
    TriMesh::Tuple mirror_edge = face_attrs[t.fid(*this)].mirror_edges[t.local_eid(*this)].value();
    assert(is_seam_edge(mirror_edge));
    TriMesh::Tuple primary_edge =
        face_attrs[mirror_edge.fid(*this)].mirror_edges[mirror_edge.local_eid(*this)].value();
    if (primary_edge.is_ccw(*this) == t.is_ccw(*this)) {
        return mirror_edge;
    } else {
        return mirror_edge.switch_vertex(*this);
    }
}

// given a seam edge with vid v retrieve the correpsonding vertex on the mirror edge
TriMesh::Tuple AdaptiveTessellation::get_mirror_vertex(const TriMesh::Tuple& t) const
{
    assert(is_seam_edge(t));
    TriMesh::Tuple mirror_edge = get_oriented_mirror_edge(t);
    return mirror_edge.switch_vertex(*this);
}

// return a vector of mirror vertices. store v itself at index 0 of the returned vector
// assume no operation has made fixed vertices outdated
std::vector<TriMesh::Tuple> AdaptiveTessellation::get_all_mirror_vertices(const TriMesh::Tuple& v)
{
    assert(is_seam_vertex(v));
    std::vector<TriMesh::Tuple> ret_vertices;
    // always put v itself at index 0
    ret_vertices.emplace_back(v);
    if (vertex_attrs[v.vid(*this)].fixed) {
        // use the same color vids to initiate Tuples
        auto same_color_uv_indices = color_to_uv_indices[uv_index_to_color[v.vid(*this)]];
        for (auto mirror_vid : same_color_uv_indices)
            ret_vertices.emplace_back(tuple_from_vertex(mirror_vid));
    } else
        ret_vertices.emplace_back(get_mirror_vertex(v));
    return ret_vertices;
}

// get all mirror_vids using navigation
std::vector<size_t> AdaptiveTessellation::get_all_mirror_vids(const TriMesh::Tuple& v)
{
    std::vector<size_t> ret_vertices_vid;
    std::queue<TriMesh::Tuple> queue;

    ret_vertices_vid.emplace_back(v.vid(*this));

    for (auto& e : get_one_ring_edges_for_vertex(v)) queue.push(e);
    while (!queue.empty()) {
        auto e = queue.front();
        queue.pop();
        if (is_seam_edge(e)) {
            auto mirror_v = get_mirror_vertex(e.switch_vertex(*this));
            if (std::find(ret_vertices_vid.begin(), ret_vertices_vid.end(), mirror_v.vid(*this)) ==
                ret_vertices_vid.end()) {
                ret_vertices_vid.emplace_back(mirror_v.vid(*this));
                for (auto& new_e : get_one_ring_edges_for_vertex(mirror_v)) queue.push(new_e);
            }
        }
    }
    return ret_vertices_vid;
}

void AdaptiveTessellation::assign_edge_curveid()
{
    for (const auto& e : get_edges()) {
        assert(e.is_valid(*this));
        if (!is_boundary_edge(e)) {
            // is an interior edge, skip ....
            continue;
        }
        // find the mid-point uv of the edge
        auto midpoint_uv =
            (vertex_attrs[e.vid(*this)].pos + vertex_attrs[e.switch_vertex(*this).vid(*this)].pos) /
            2.;
        // use the mid-point uv to find edge curve id
        int curve_id = -1;
        double t = 0.;
        std::tie(curve_id, t) = mesh_parameters.m_boundary.uv_to_t(midpoint_uv);
        // assign the curve id to the edge
        edge_attrs[e.eid(*this)].curve_id = std::make_optional<int>(curve_id);
    }
}

} // namespace adaptive_tessellation
