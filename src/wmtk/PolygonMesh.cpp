#include "PolygonMesh.hpp"

namespace wmtk {
PolygonMesh::PolygonMesh()
    : Mesh(2)
    , m_next_handle(register_attribute<long>("m_next", PrimitiveType::HalfEdge, 1))
    , m_prev_handle(register_attribute<long>("m_prev", PrimitiveType::HalfEdge, 1))
    , m_to_handle(register_attribute<long>("m_to", PrimitiveType::HalfEdge, 1))
    , m_out_handle(register_attribute<long>("m_out", PrimitiveType::Vertex, 1))
    , m_hf_handle(register_attribute<long>("m_hf", PrimitiveType::HalfEdge, 1))
    , m_fh_handle(register_attribute<long>("m_fh", PrimitiveType::Face, 1))
    , m_f_is_hole_handle(register_attribute<char>("m_f_is_hole", PrimitiveType::Face, 1))
{}
PolygonMesh::PolygonMesh(const PolygonMesh& o) = default;
PolygonMesh::PolygonMesh(PolygonMesh&& o) = default;
PolygonMesh& PolygonMesh::operator=(const PolygonMesh& o) = default;
PolygonMesh& PolygonMesh::operator=(PolygonMesh&& o) = default;

Tuple PolygonMesh::split_edge(const Tuple& t, Accessor<long>& hash_accessor)
{
    // TODO Need to implement
    assert(false);
    return t;
}

Tuple PolygonMesh::collapse_edge(const Tuple& t, Accessor<long>& hash_accessor)
{
    // TODO Need to implement
    assert(false);
    return t;
}

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
        long edge_hid = h % 2; // Local id for the halfedge with implicit opp

        return Tuple(
            new_local_vid,
            tuple.m_local_eid,
            tuple.m_local_fid,
            (edge_hid == 0) ? (h + 1) : (h - 1),
            tuple.m_hash);
    }
    case PrimitiveType::Tetrahedron:
    default: throw std::runtime_error("Tuple switch: Invalid primitive type"); break;
    }
}

Tuple PolygonMesh::next_halfedge(const Tuple& h_tuple) const
{
    ConstAccessor<long> next_accessor = create_const_accessor<long>(m_next_handle);
    long next_h = next_accessor.scalar_attribute(h_tuple);

    return Tuple(
        h_tuple.m_local_vid,
        h_tuple.m_local_eid,
        h_tuple.m_local_fid,
        next_h,
        h_tuple.m_hash);
}

Tuple PolygonMesh::prev_halfedge(const Tuple& h_tuple) const
{
    ConstAccessor<long> prev_accessor = create_const_accessor<long>(m_prev_handle);
    long prev_h = prev_accessor.scalar_attribute(h_tuple);

    return Tuple(
        h_tuple.m_local_vid,
        h_tuple.m_local_eid,
        h_tuple.m_local_fid,
        prev_h,
        h_tuple.m_hash);
}

Tuple PolygonMesh::opp_halfedge(const Tuple& h_tuple) const
{
    long h = h_tuple.m_global_cid;
    long edge_hid = h % 2; // Local id for the halfedge with implicit opp
    return Tuple(
        h_tuple.m_local_vid,
        h_tuple.m_local_eid,
        h_tuple.m_local_fid,
        (edge_hid == 0) ? (h + 1) : (h - 1),
        h_tuple.m_hash);
}

bool PolygonMesh::is_connectivity_valid() const
{
    // TODO
    assert(false);
    return true;
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
    case PrimitiveType::Tetrahedron:
    default: throw std::runtime_error("Invalid primitive type"); break;
    }
}

bool PolygonMesh::is_ccw(const Tuple& tuple) const
{
    return (tuple.m_local_vid == 0);
}

bool PolygonMesh::is_hole_face(const Tuple& f_tuple) const
{
    ConstAccessor<char> f_is_hole_accessor = create_const_accessor<char>(m_f_is_hole_handle);
    auto f_is_hole = f_is_hole_accessor.scalar_attribute(f_tuple);
    return (f_is_hole == 1);
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
        v_tuple = switch_edge(switch_face(v_tuple));
    } while (v_tuple != tuple);

    return false;
}

bool PolygonMesh::is_boundary_edge(const Tuple& tuple) const
{
    return (is_hole_face(tuple)) || (is_hole_face(opp_halfedge(tuple)));
}

bool PolygonMesh::is_valid(const Tuple& tuple, ConstAccessor<long>& hash_accessor) const
{
    // TODO
    assert(false);
    return true;
}

void PolygonMesh::initialize(
    Eigen::Ref<const VectorXl> next,
    Eigen::Ref<const VectorXl> boundary_faces)
{
    // TODO
    assert(false);
}

void PolygonMesh::initialize(Eigen::Ref<const RowVectors3l> F)
{
    // TODO
    assert(false);
}

long PolygonMesh::id(const Tuple& tuple, PrimitiveType type) const
{
    switch (type) {
    case PrimitiveType::Vertex: {
        ConstAccessor<long> to_accessor = create_const_accessor<long>(m_to_handle);
        if (tuple.m_local_vid == 0) {
            return to_accessor.scalar_attribute(tuple);
        } else {
            return to_accessor.scalar_attribute(opp_halfedge(tuple));
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
    case PrimitiveType::Tetrahedron:
    default: throw std::runtime_error("Invalid primitive type"); break;
    }
}

} // namespace wmtk