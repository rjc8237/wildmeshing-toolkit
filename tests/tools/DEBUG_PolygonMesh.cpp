#include "DEBUG_PolygonMesh.hpp"
#include <queue>
#include <stdexcept>
#include <wmtk/utils/Logger.hpp>

namespace wmtk::tests {

namespace {
// Implicit opposite map for halfedges paired as [e] = {2*e, 2*e + 1}
long implicit_opp(long h)
{
    return ((h % 2) == 0) ? (h + 1) : (h - 1);
}
} // namespace

DEBUG_PolygonMesh::DEBUG_PolygonMesh(const PolygonMesh& m)
    : PolygonMesh(m)
{}
DEBUG_PolygonMesh::DEBUG_PolygonMesh(PolygonMesh&& m)
    : PolygonMesh(std::move(m))
{}


bool DEBUG_PolygonMesh::operator==(const DEBUG_PolygonMesh& o) const
{
    return static_cast<const PolygonMesh&>(*this) == static_cast<const PolygonMesh&>(o);
}
bool DEBUG_PolygonMesh::operator!=(const DEBUG_PolygonMesh& o) const
{
    return !(*this == o);
}

long DEBUG_PolygonMesh::id(const Tuple& tuple, PrimitiveType type) const
{
    return PolygonMesh::id(tuple, type);
}

const MeshAttributeHandle<long>& DEBUG_PolygonMesh::next_handle() const
{
    return m_next_handle;
}

Tuple DEBUG_PolygonMesh::halfedge_tuple_from_vertex_in_face(long vid, long fid) const
{
    Tuple f_tuple = tuple_from_id(PrimitiveType::Face, fid);
    Tuple f_tuple_iter = f_tuple;
    do {
        if (id(f_tuple_iter, PrimitiveType::Vertex) == vid) {
            break;
        }
        f_tuple_iter = next_halfedge(f_tuple_iter);
    } while (id(f_tuple_iter, PrimitiveType::HalfEdge) != id(f_tuple, PrimitiveType::HalfEdge));

    assert(id(f_tuple_iter, PrimitiveType::Vertex) == vid);
    return f_tuple_iter;
}

std::array<long, 4> DEBUG_PolygonMesh::primitive_counts() const
{
    return std::array<long, 4>{
        {long(get_all(PrimitiveType::Vertex).size()),
         long(get_all(PrimitiveType::Edge).size()),
         long(get_all(PrimitiveType::Face).size()),
         long(get_all(PrimitiveType::HalfEdge).size())}};
}

long DEBUG_PolygonMesh::find_simply_connected_components(std::vector<long>& component_ids) const
{
    component_ids.resize(capacity(PrimitiveType::HalfEdge), -1);
    const auto all_halfedge_tuples = get_all(PrimitiveType::HalfEdge);
    long component_id = 0;

    // BFS
    for (auto halfedge_tuple : all_halfedge_tuples) {
        const long hid = id(halfedge_tuple, PrimitiveType::HalfEdge);
        if (component_ids[hid] != -1) {
            continue;
        } // visited

        std::queue<long> q;
        q.push(hid);
        while (!q.empty()) {
            long cur_hid = q.front();
            q.pop();

            component_ids[cur_hid] = component_id;
            auto h_tuple = tuple_from_id(PrimitiveType::HalfEdge, cur_hid);

            // Add opposite halfedge to component
            auto opp_tuple = opp_halfedge(h_tuple);
            long opp_hid = id(opp_tuple, PrimitiveType::HalfEdge);
            if (component_ids[opp_hid] == -1) {
                q.push(opp_hid);
            }

            // Add next halfedge to component
            auto next_tuple = next_halfedge(h_tuple);
            long next_hid = id(next_tuple, PrimitiveType::HalfEdge);
            if (component_ids[next_hid] == -1) {
                q.push(next_hid);
            }
        }

        component_id++;
    }

    return component_id;
}

long DEBUG_PolygonMesh::count_simply_connected_components() const
{
    std::vector<long> component_ids;
    return find_simply_connected_components(component_ids);
}

long DEBUG_PolygonMesh::count_boundary_loops() const
{
    std::vector<long> boundary_loop_ids(capacity(PrimitiveType::Edge), -1);
    const auto all_edge_tuples = get_all(PrimitiveType::Edge);
    long boundary_loop_id = 0;
    for (auto edge_tuple : all_edge_tuples) {
        if (!is_boundary(edge_tuple)) {
            continue;
        }
        const long eid = id(edge_tuple, PrimitiveType::Edge);
        if (boundary_loop_ids[eid] != -1) {
            continue;
        } // visited

        Tuple cur_edge = edge_tuple;
        long cur_eid = id(cur_edge, PrimitiveType::Edge);
        // find one boundary loop
        while (boundary_loop_ids[eid] == -1 || cur_eid != eid) {
            boundary_loop_ids[cur_eid] = boundary_loop_id;

            // find next boundary edge
            cur_edge = switch_edge(switch_vertex(cur_edge));
            while (!is_boundary(cur_edge)) {
                cur_edge = switch_edge(switch_face(cur_edge));
            }

            cur_eid = id(cur_edge, PrimitiveType::Edge);
        }

        boundary_loop_id++;
    }

    return boundary_loop_id;
}

long DEBUG_PolygonMesh::count_hole_faces() const
{
    auto faces = get_all(PrimitiveType::Face);
    long num_hole_faces = 0;
    for (const auto& face : faces) {
        if (is_hole_face(face)) {
            ++num_hole_faces;
        }
    }

    return num_hole_faces;
}

long DEBUG_PolygonMesh::euler_characteristic() const
{
    auto counts = primitive_counts();
    long num_boundary_loops = count_boundary_loops();
    long num_hole_faces = count_hole_faces();
    return counts[0] - counts[1] + counts[2] + num_boundary_loops - num_hole_faces;
}

long DEBUG_PolygonMesh::genus() const
{
    if (count_simply_connected_components() != 1) {
        wmtk::logger().error("not a connected surface\n");
        throw std::runtime_error("GenusComputeError");
    }
    return (2 - euler_characteristic()) / 2;
}

bool DEBUG_PolygonMesh::is_connectivity_valid() const
{
    ConstAccessor<char> v_flag_accessor = get_flag_accessor(PrimitiveType::Vertex);
    ConstAccessor<char> e_flag_accessor = get_flag_accessor(PrimitiveType::Edge);
    ConstAccessor<char> f_flag_accessor = get_flag_accessor(PrimitiveType::Face);
    ConstAccessor<char> h_flag_accessor = get_flag_accessor(PrimitiveType::HalfEdge);

    // Halfedge connectivity
    long n_halfedges = 0;
    for (long hid = 0; hid < capacity(PrimitiveType::HalfEdge); ++hid) {
        if (get_index_access(h_flag_accessor).scalar_attribute(hid) == 0) {
            continue;
        } else {
            ++n_halfedges;
        }

        if (!is_halfedge_connectivity_valid(hid)) {
            wmtk::logger().error("halfedge connectivity is not valid for halfedge", hid);
            return false;
        }
    }

    // Vertex Connectivity
    // TODO Check vertex boundary conditions
    long n_vertices = 0;
    for (long vid = 0; vid < capacity(PrimitiveType::Vertex); ++vid) {
        // Count vertices
        if (get_index_access(v_flag_accessor).scalar_attribute(vid) == 0) {
            continue;
        } else {
            ++n_vertices;
        }

        // Check local vertex validity
        if (!is_vertex_connectivity_valid(vid)) {
            wmtk::logger().error("vertex connectivity is not valid for vertex {}", vid);
            return false;
        }
    }

    // Edge Connectivity
    long n_edges = 0;
    for (long eid = 0; eid < capacity(PrimitiveType::Edge); ++eid) {
        // Count edges
        if (get_index_access(e_flag_accessor).scalar_attribute(eid) == 0) {
            continue;
        } else {
            ++n_edges;
        }

        // Check local edge validity
        if (!is_edge_connectivity_valid(eid)) {
            wmtk::logger().error("edge connectivity is not valid for edge {}", eid);
            return false;
        }
    }

    // Face Connectivity
    long n_faces = 0;
    for (long fid = 0; fid < capacity(PrimitiveType::Face); ++fid) {
        // Count faces
        if (get_index_access(f_flag_accessor).scalar_attribute(fid) == 0) {
            continue;
        } else {
            ++n_faces;
        }

        // Check local face validity
        if (!is_face_connectivity_valid(fid)) {
            wmtk::logger().error("face connectivity is not valid for face {}", fid);
            return false;
        }
    }

    // Check element counts are consistent
    if (n_halfedges != (2 * n_edges)) {
        wmtk::logger().error(
            "edge count {} and halfedge count {} are inconsistent",
            n_edges,
            n_halfedges);
        return false;
    }
    long n_circulator_orbits = count_vertices_from_orbits();
    if (n_vertices != n_circulator_orbits) {
        wmtk::logger().error(
            "vertex count {} and circulator orbit count {} are inconsistent",
            n_vertices,
            n_circulator_orbits);
        return false;
    }
    long n_next_orbits = count_faces_from_orbits();
    if (n_faces != n_next_orbits) {
        wmtk::logger().error(
            "face count {} and next orbit count {} are inconsistent",
            n_faces,
            n_next_orbits);
        return false;
    }

    return true;
}

bool DEBUG_PolygonMesh::is_local_connectivity_valid(const Tuple& tuple, PrimitiveType type) const
{
    switch (type) {
    case PrimitiveType::Vertex: {
        return is_vertex_connectivity_valid(id(tuple, type));
    }
    case PrimitiveType::Edge: {
        return is_edge_connectivity_valid(id(tuple, type));
    }
    case PrimitiveType::Face: {
        return is_face_connectivity_valid(id(tuple, type));
    }
    case PrimitiveType::HalfEdge: {
        return is_halfedge_connectivity_valid(id(tuple, type));
    }
    case PrimitiveType::Tetrahedron:
    default: throw std::runtime_error("Invalid primitive type"); break;
    }
}

bool DEBUG_PolygonMesh::is_vertex_connectivity_valid(long vid) const
{
    // Vertex id is in range
    long v_cap = capacity(PrimitiveType::Vertex);
    if ((vid < 0) || (vid >= v_cap)) {
        wmtk::logger().error("vertex id {} is out of range [0, {})", vid, v_cap);
        return false;
    }

    // Vertex is not deleted
    ConstAccessor<char> v_flag_accessor = get_flag_accessor(PrimitiveType::Vertex);
    if (get_index_access(v_flag_accessor).scalar_attribute(vid) == 0) {
        wmtk::logger().error("vertex {} is deleted", vid);
        return false;
    }

    // The circulator orbit of the vertex halfedge has constant vertex index fid
    ConstAccessor<long> out_accessor = create_const_accessor<long>(m_out_handle);
    ConstAccessor<long> to_accessor = create_const_accessor<long>(m_to_handle);
    ConstAccessor<long> prev_accessor = create_const_accessor<long>(m_prev_handle);
    long hid_start = implicit_opp(get_index_access(out_accessor).scalar_attribute(vid));
    long hid_iter = hid_start;
    long n_halfedges = capacity(PrimitiveType::HalfEdge);
    long n_circulations = 0;
    do {
        long hid_to = get_index_access(to_accessor).scalar_attribute(hid_iter);
        if (hid_to != vid) {
            wmtk::logger().error(
                "halfedge {} in the circulator orbit of vertex {} points to {}",
                hid_iter,
                vid,
                hid_to);
            return false;
        }
        if (n_circulations > n_halfedges) {
            wmtk::logger().error("more vertex circulations than possible attempted");
            return false;
        }

        // Circulate halfedge around vertex
        hid_iter = get_index_access(prev_accessor).scalar_attribute(implicit_opp(hid_iter));
        ++n_circulations;
    } while (hid_iter != hid_start);

    return true;
}

bool DEBUG_PolygonMesh::is_edge_connectivity_valid(long eid) const
{
    // Edge id is in range
    long e_cap = capacity(PrimitiveType::Edge);
    if ((eid < 0) || (eid >= e_cap)) {
        wmtk::logger().error("edge id {} is out of range [0, {})", eid, e_cap);
        return false;
    }

    // Edge is not deleted
    ConstAccessor<char> e_flag_accessor = get_flag_accessor(PrimitiveType::Edge);
    if (get_index_access(e_flag_accessor).scalar_attribute(eid) == 0) {
        wmtk::logger().error("edge {} is deleted", eid);
        return false;
    }

    // opp and the edge to halfedge and halfedge to edge maps are implicit
    // We simply check that the corresponding halfedge ids 2e and 2e + 1 are in range
    long h_cap = capacity(PrimitiveType::HalfEdge);
    if (2 * eid > h_cap) {
        wmtk::logger().error("edge {} is inconsistent with halfedge range [0, {})", eid, h_cap);
        return false;
    }

    return true;
}

bool DEBUG_PolygonMesh::is_face_connectivity_valid(long fid) const
{
    // Face id is in range
    long f_cap = capacity(PrimitiveType::Face);
    if ((fid < 0) || (fid >= f_cap)) {
        wmtk::logger().error("face id {} is out of range [0, {})", fid, f_cap);
        return false;
    }

    // Face is not deleted
    ConstAccessor<char> f_flag_accessor = get_flag_accessor(PrimitiveType::Face);
    if (get_index_access(f_flag_accessor).scalar_attribute(fid) == 0) {
        wmtk::logger().error("face {} is deleted", fid);
        return false;
    }

    // The next orbit of the face halfedge has constant face index fid
    ConstAccessor<long> fh_accessor = create_const_accessor<long>(m_fh_handle);
    ConstAccessor<long> hf_accessor = create_const_accessor<long>(m_hf_handle);
    ConstAccessor<long> next_accessor = create_const_accessor<long>(m_next_handle);
    long hid_start = get_index_access(fh_accessor).scalar_attribute(fid);
    long hid_iter = hid_start;
    long n_halfedges = capacity(PrimitiveType::HalfEdge);
    long n_circulations = 0;
    do {
        long hid_face = get_index_access(hf_accessor).scalar_attribute(hid_iter);
        if (hid_face != fid) {
            wmtk::logger().error(
                "halfedge {} in the next orbit of face {} is adjacent to {}",
                hid_iter,
                fid,
                hid_face);
            return false;
        }
        if (n_circulations > n_halfedges) {
            wmtk::logger().error("more face circulations than possible attempted");
            return false;
        }

        // Circulate halfedge around face
        hid_iter = get_index_access(next_accessor).scalar_attribute(hid_iter);
        ++n_circulations;
    } while (hid_iter != hid_start);

    return true;
}

bool DEBUG_PolygonMesh::is_halfedge_connectivity_valid(long hid) const
{
    // Halfedge id is in range
    long h_cap = capacity(PrimitiveType::HalfEdge);
    if ((hid < 0) || (hid >= h_cap)) {
        wmtk::logger().error("halfedge {} is out of range [0, {})", hid, h_cap);
        return false;
    }

    // Halfedge is not deleted
    ConstAccessor<char> h_flag_accessor = get_flag_accessor(PrimitiveType::HalfEdge);
    if (get_index_access(h_flag_accessor).scalar_attribute(hid) == 0) {
        wmtk::logger().error("halfedge {} is deleted", hid);
        return false;
    }

    // next and prev are (locally) inverse
    ConstAccessor<long> next_accessor = create_const_accessor<long>(m_next_handle);
    ConstAccessor<long> prev_accessor = create_const_accessor<long>(m_prev_handle);
    long hid_next = get_index_access(next_accessor).scalar_attribute(hid);
    long hid_np = get_index_access(prev_accessor).scalar_attribute(hid_next);
    if (hid_np != hid) {
        wmtk::logger().error("halfedge {} is next of {} but has prev {}", hid_next, hid, hid_np);
        return false;
    }
    long hid_prev = get_index_access(prev_accessor).scalar_attribute(hid);
    long hid_pn = get_index_access(next_accessor).scalar_attribute(hid_prev);
    if (hid_pn != hid) {
        wmtk::logger().error("halfedge {} is prev of {} but has next {}", hid_prev, hid, hid_pn);
        return false;
    }

    return true;
}

long DEBUG_PolygonMesh::count_vertices_from_orbits() const
{
    ConstAccessor<char> h_flag_accessor = get_flag_accessor(PrimitiveType::HalfEdge);
    ConstAccessor<long> prev_accessor = create_const_accessor<long>(m_prev_handle);

    long n_halfedges = capacity(PrimitiveType::HalfEdge);
    std::vector<bool> visited(n_halfedges, false);
    std::vector<std::vector<long>> vertex_cycles(0);
    vertex_cycles.reserve(n_halfedges);
    for (long hid = 0; hid < n_halfedges; ++hid) {
        // Skip deleted halfedges
        if (get_index_access(h_flag_accessor).scalar_attribute(hid) == 0) {
            continue;
        }

        // Build the vertex orbit of the halfedge under the circulator if it hasn't been seen yet
        if (!visited[hid]) {
            vertex_cycles.emplace_back(std::vector<long>());
            long hid_iter = hid;
            while (true) {
                visited[hid_iter] = true;
                vertex_cycles.back().push_back(hid_iter);
                hid_iter = implicit_opp(hid_iter);
                hid_iter = get_index_access(prev_accessor).scalar_attribute(hid_iter);
                if (hid_iter == hid) {
                    break;
                }
            }
        }
    }

    return vertex_cycles.size();
}

long DEBUG_PolygonMesh::count_faces_from_orbits() const
{
    ConstAccessor<char> h_flag_accessor = get_flag_accessor(PrimitiveType::HalfEdge);
    ConstAccessor<long> next_accessor = create_const_accessor<long>(m_next_handle);

    long n_halfedges = capacity(PrimitiveType::HalfEdge);
    std::vector<bool> visited(n_halfedges, false);
    std::vector<std::vector<long>> face_cycles(0);
    face_cycles.reserve(n_halfedges);
    for (long hid = 0; hid < n_halfedges; ++hid) {
        // Skip deleted halfedges
        if (get_index_access(h_flag_accessor).scalar_attribute(hid) == 0) {
            continue;
        }

        // Build the face orbit of the halfedge under next if it hasn't been seen yet
        if (!visited[hid]) {
            face_cycles.emplace_back(std::vector<long>());
            long hid_iter = hid;
            while (true) {
                visited[hid_iter] = true;
                face_cycles.back().push_back(hid_iter);
                hid_iter = get_index_access(next_accessor).scalar_attribute(hid_iter);
                if (hid_iter == hid) {
                    break;
                }
            }
        }
    }

    return face_cycles.size();
}

} // namespace wmtk::tests
