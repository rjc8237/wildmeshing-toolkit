#include "PolygonMesh.hpp"
#include <spdlog/spdlog.h>
#include <wmtk/utils/polygon_mesh_topology_initialization.h>
#include <wmtk/utils/Logger.hpp>
#include <wmtk/utils/polygon_mesh_topology_validity.hpp>

namespace wmtk {

namespace {
// Implicit opposite map for halfedges paired as [e] = {2*e, 2*e + 1}
long implicit_opp(long h)
{
    return ((h % 2) == 0) ? (h + 1) : (h - 1);
}
} // namespace

PolygonMesh::PolygonMesh()
    : Mesh(2, 4, PrimitiveType::HalfEdge)
    , m_next_handle(register_attribute<long>("m_next", PrimitiveType::HalfEdge, 1))
    , m_prev_handle(register_attribute<long>("m_prev", PrimitiveType::HalfEdge, 1))
    , m_to_handle(register_attribute<long>("m_to", PrimitiveType::HalfEdge, 1))
    , m_out_handle(register_attribute<long>("m_out", PrimitiveType::Vertex, 1))
    , m_hf_handle(register_attribute<long>("m_hf", PrimitiveType::HalfEdge, 1))
    , m_fh_handle(register_attribute<long>("m_fh", PrimitiveType::Face, 1))
    , m_f_is_hole_handle(register_attribute<char>("m_f_is_hole", PrimitiveType::Face, 1))
{
    assert(get_max_primitive_type_id({PrimitiveType::Face, PrimitiveType::HalfEdge}) == 4);
}

PolygonMesh::PolygonMesh(const PolygonMesh& o) = default;
PolygonMesh::PolygonMesh(PolygonMesh&& o) = default;
PolygonMesh& PolygonMesh::operator=(const PolygonMesh& o) = default;
PolygonMesh& PolygonMesh::operator=(PolygonMesh&& o) = default;

Tuple PolygonMesh::switch_tuple(const Tuple& tuple, PrimitiveType type) const
{
    long new_local_vid = (tuple.m_local_vid + 1) % 2; // All switches swap local vid

    switch (type) {
    case PrimitiveType::Vertex:
        return Tuple(
            new_local_vid,
            tuple.m_local_eid,
            tuple.m_local_fid,
            tuple.m_global_cid,
            tuple.m_hash);

    case PrimitiveType::Edge: {
        long h = tuple.m_global_cid;

        // Switched edge is prev for local vertex index 0 and next for local index 1
        auto e_switch_handle = (tuple.m_local_vid == 0) ? m_prev_handle : m_next_handle;
        ConstAccessor<long> e_switch_accessor = create_const_accessor<long>(e_switch_handle);
        long e_switch = e_switch_accessor.scalar_attribute(tuple);

        return Tuple(new_local_vid, tuple.m_local_eid, tuple.m_local_fid, e_switch, tuple.m_hash);
    }
    case PrimitiveType::Face: {
        long h = tuple.m_global_cid;

        return Tuple(
            new_local_vid,
            tuple.m_local_eid,
            tuple.m_local_fid,
            implicit_opp(h),
            tuple.m_hash);
    }
    case PrimitiveType::HalfEdge:
    case PrimitiveType::Tetrahedron:
    default: throw std::runtime_error("Tuple switch: Invalid primitive type"); break;
    }
}

Tuple PolygonMesh::tuple_from_id(const PrimitiveType type, const long gid) const
{
    // TODO Need to get hashes

    switch (type) {
    case PrimitiveType::Vertex: {
        ConstAccessor<long> out_accessor = create_const_accessor<long>(m_out_handle);
        auto h = out_accessor.index_access().scalar_attribute(gid);
        return Tuple(0, -1, -1, h, 0);
    }
    case PrimitiveType::Edge: {
        long h = 2 * gid;
        return Tuple(0, -1, -1, h, 0);
    }
    case PrimitiveType::Face: {
        ConstAccessor<long> fh_accessor = create_const_accessor<long>(m_fh_handle);
        auto h = fh_accessor.index_access().scalar_attribute(gid);
        return Tuple(0, -1, -1, h, 0);
    }
    case PrimitiveType::HalfEdge: {
        return Tuple(0, -1, -1, gid, 0);
    }
    case PrimitiveType::Tetrahedron:
    default: throw std::runtime_error("Invalid primitive type"); break;
    }
}

bool PolygonMesh::is_connectivity_valid() const
{
    // Validity is ensured by construction
    return true;
}

bool PolygonMesh::is_valid(const Tuple& tuple, ConstAccessor<long>& hash_accessor) const
{
    if (tuple.is_null()) return false;
    if (tuple.m_global_cid < 0) {
        return false;
    }

    // Check if local vid is in {0, 1}
    if ((tuple.m_local_vid != 0) && (tuple.m_local_vid != 1)) {
        return false;
    }

    return Mesh::is_hash_valid(tuple, hash_accessor);
}

bool PolygonMesh::is_ccw(const Tuple& tuple) const
{
    return (tuple.m_local_vid == 0);
}

bool PolygonMesh::is_boundary(const Tuple& tuple) const
{
    return is_boundary_edge(tuple);
}

bool PolygonMesh::is_boundary_vertex(const Tuple& tuple) const
{
    Tuple v_tuple = tuple;
    do {
        if (is_boundary(v_tuple)) {
            return true;
        }
        v_tuple = switch_tuples(v_tuple, {PrimitiveType::Face, PrimitiveType::Edge});
    } while (v_tuple != tuple);

    return false;
}

bool PolygonMesh::is_boundary_edge(const Tuple& tuple) const
{
    return (is_hole_face(tuple)) || (is_hole_face(opp_halfedge(tuple)));
}

Tuple PolygonMesh::next_halfedge(const Tuple& h_tuple) const
{
    return switch_tuples(h_tuple, {PrimitiveType::Vertex, PrimitiveType::Edge});
}

Tuple PolygonMesh::prev_halfedge(const Tuple& h_tuple) const
{
    return switch_tuples(h_tuple, {PrimitiveType::Edge, PrimitiveType::Vertex});
}

Tuple PolygonMesh::opp_halfedge(const Tuple& h_tuple) const
{
    long h = h_tuple.m_global_cid;
    return Tuple(
        h_tuple.m_local_vid,
        h_tuple.m_local_eid,
        h_tuple.m_local_fid,
        implicit_opp(h),
        h_tuple.m_hash);
}

bool PolygonMesh::is_hole_face(const Tuple& f_tuple) const
{
    ConstAccessor<char> f_is_hole_accessor = create_const_accessor<char>(m_f_is_hole_handle);
    auto f_is_hole = f_is_hole_accessor.scalar_attribute(f_tuple);
    return (f_is_hole == 1);
}

void PolygonMesh::initialize(
    Eigen::Ref<const VectorXl> next,
    Eigen::Ref<const VectorXl> prev,
    Eigen::Ref<const VectorXl> to,
    Eigen::Ref<const VectorXl> out,
    Eigen::Ref<const VectorXl> he2f,
    Eigen::Ref<const VectorXl> f2he)
{
    assert(utils::are_polygon_mesh_edges_valid(next, prev));
    assert(utils::are_polygon_mesh_vertices_valid(prev, to, out));
    assert(utils::are_polygon_mesh_faces_valid(next, he2f, f2he));

    std::vector<long> cap{
        static_cast<long>(out.size()),
        static_cast<long>(next.size() / 2),
        static_cast<long>(f2he.rows()),
        0,
        static_cast<long>(next.size())};
    set_capacities(cap);

    // get Accessors for topology
    Accessor<long> next_accessor = create_accessor<long>(m_next_handle);
    Accessor<long> prev_accessor = create_accessor<long>(m_prev_handle);
    Accessor<long> to_accessor = create_accessor<long>(m_to_handle);
    Accessor<long> out_accessor = create_accessor<long>(m_out_handle);
    Accessor<long> hf_accessor = create_accessor<long>(m_hf_handle);
    Accessor<long> fh_accessor = create_accessor<long>(m_fh_handle);
    Accessor<char> f_is_hole_accessor = create_accessor<char>(m_f_is_hole_handle);

    Accessor<char> v_flag_accessor = get_flag_accessor(PrimitiveType::Vertex);
    Accessor<char> e_flag_accessor = get_flag_accessor(PrimitiveType::Edge);
    Accessor<char> f_flag_accessor = get_flag_accessor(PrimitiveType::Face);
    Accessor<char> h_flag_accessor = get_flag_accessor(PrimitiveType::HalfEdge);

    // iterate over the vectors and fill attributes
    for (long vi = 0; vi < capacity(PrimitiveType::Vertex); ++vi) {
        v_flag_accessor.index_access().scalar_attribute(vi) |= 0x1;
        out_accessor.index_access().scalar_attribute(vi) = out(vi);
    }
    for (long ei = 0; ei < capacity(PrimitiveType::Edge); ++ei) {
        e_flag_accessor.index_access().scalar_attribute(ei) |= 0x1;
    }
    for (long fi = 0; fi < capacity(PrimitiveType::Face); ++fi) {
        f_flag_accessor.index_access().scalar_attribute(fi) |= 0x1;
        fh_accessor.index_access().scalar_attribute(fi) = f2he(fi);
        f_is_hole_accessor.index_access().scalar_attribute(fi) = 0;
    }
    for (long hi = 0; hi < capacity(PrimitiveType::HalfEdge); ++hi) {
        h_flag_accessor.index_access().scalar_attribute(hi) |= 0x1;
        next_accessor.index_access().scalar_attribute(hi) = next(hi);
        prev_accessor.index_access().scalar_attribute(hi) = prev(hi);
        to_accessor.index_access().scalar_attribute(hi) = to(hi);
        hf_accessor.index_access().scalar_attribute(hi) = he2f(hi);
    }
}

void PolygonMesh::initialize(Eigen::Ref<const VectorXl> next)
{
    auto [prev, to, out, he2f, f2he] = utils::polygon_mesh_topology_initialization(next);
    initialize(next, prev, to, out, he2f, f2he);
}

void PolygonMesh::initialize_fv(const std::vector<std::vector<long>>& F)
{
    auto [next, prev, to, out, he2f, f2he, bnd_loops] =
        utils::polygon_mesh_fv_topology_initialization(F);
    initialize(next, prev, to, out, he2f, f2he);

    Accessor<char> f_is_hole_accessor = create_accessor<char>(m_f_is_hole_handle);
    for (const auto& bnd_loop : bnd_loops) {
        f_is_hole_accessor.index_access().scalar_attribute(he2f(bnd_loop)) = 1;
    }
}

void PolygonMesh::initialize_fv(Eigen::Ref<const RowVectors3l> F)
{
    auto [next, prev, to, out, he2f, f2he, bnd_loops] =
        utils::polygon_mesh_fv_topology_initialization(F);
    initialize(next, prev, to, out, he2f, f2he);

    Accessor<char> f_is_hole_accessor = create_accessor<char>(m_f_is_hole_handle);
    for (const auto& bnd_loop : bnd_loops) {
        f_is_hole_accessor.index_access().scalar_attribute(he2f(bnd_loop)) = 1;
    }
}


long PolygonMesh::id(const Tuple& tuple, PrimitiveType type) const
{
    switch (type) {
    case PrimitiveType::Vertex: {
        ConstAccessor<long> to_accessor = create_const_accessor<long>(m_to_handle);
        if (tuple.m_local_vid == 0) {
            return to_accessor.scalar_attribute(opp_halfedge(tuple));
        } else {
            return to_accessor.scalar_attribute(tuple);
        }
    }
    case PrimitiveType::Edge: {
        long h = tuple.m_global_cid;
        return h / 2;
    }
    case PrimitiveType::Face: {
        ConstAccessor<long> hf_accessor = create_const_accessor<long>(m_hf_handle);
        long f = hf_accessor.scalar_attribute(tuple);
        return f;
    }
    case PrimitiveType::HalfEdge: {
        return tuple.m_global_cid;
    }
    case PrimitiveType::Tetrahedron:
    default: throw std::runtime_error("Invalid primitive type"); break;
    }
}

} // namespace wmtk