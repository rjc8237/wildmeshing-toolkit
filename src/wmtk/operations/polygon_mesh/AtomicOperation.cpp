#include "AtomicOperation.hpp"

namespace wmtk::operations::polygon_mesh {

namespace {
// Implicit opposite map for halfedges paired as [e] = {2*e, 2*e + 1}
long implicit_opp(long h)
{
    return ((h % 2) == 0) ? (h + 1) : (h - 1);
}
} // namespace

AtomicOperation::AtomicOperation(PolygonMesh& m)
    : PolygonMeshOperation(m)
    , m_hn_accessor(m.create_accessor<long>(m.m_hn_handle))
    , m_hp_accessor(m.create_accessor<long>(m.m_hp_handle))
    , m_hf_accessor(m.create_accessor<long>(m.m_hf_handle))
    , m_hv_accessor(m.create_accessor<long>(m.m_hv_handle))
    , m_fh_accessor(m.create_accessor<long>(m.m_fh_handle))
    , m_vh_accessor(m.create_accessor<long>(m.m_vh_handle))
    , m_v_flag_accessor(m.get_flag_accessor(PrimitiveType::Vertex))
    , m_e_flag_accessor(m.get_flag_accessor(PrimitiveType::Edge))
    , m_f_flag_accessor(m.get_flag_accessor(PrimitiveType::Face))
    , m_h_flag_accessor(m.get_flag_accessor(PrimitiveType::HalfEdge))
{}

AtomicOperation::AtomicOperation(Mesh& m)
    : AtomicOperation(dynamic_cast<PolygonMesh&>(m))
{}

long AtomicOperation::new_face(bool is_hole)
{
    mesh().reserve_attributes(PrimitiveType::Face, mesh().capacity(PrimitiveType::Face) + 1);
    long f_id = mesh().request_simplex_indices(PrimitiveType::Face, 1)[0];
    return is_hole ? -f_id : f_id;
}

long AtomicOperation::new_edge()
{
    mesh().reserve_attributes(PrimitiveType::Edge, mesh().capacity(PrimitiveType::Edge) + 1);
    mesh().reserve_attributes(PrimitiveType::HalfEdge, mesh().capacity(PrimitiveType::Edge) + 2);
    return mesh().request_simplex_indices(PrimitiveType::Edge, 1)[0];
}

long AtomicOperation::new_vertex()
{
    mesh().reserve_attributes(PrimitiveType::Vertex, mesh().capacity(PrimitiveType::Vertex) + 1);
    return mesh().request_simplex_indices(PrimitiveType::Vertex, 1)[0];
}

void AtomicOperation::delete_face(long face_id)
{
    m_f_flag_accessor.index_access().scalar_attribute(face_id) = 0;
}

void AtomicOperation::delete_edge(long edge_id)
{
    m_e_flag_accessor.index_access().scalar_attribute(edge_id) = 0;
    m_h_flag_accessor.index_access().scalar_attribute(2 * edge_id) = 0;
    m_h_flag_accessor.index_access().scalar_attribute(2 * edge_id + 1) = 0;
}

void AtomicOperation::delete_vertex(long vertex_id)
{
    m_v_flag_accessor.index_access().scalar_attribute(vertex_id) = 0;
}

void AtomicOperation::set_face(long halfedge_id, long face_id)
{
    long h_iter_id = halfedge_id;
    do {
        m_hf_accessor.index_access().scalar_attribute(h_iter_id) = face_id;
        h_iter_id = m_hn_accessor.index_access().scalar_attribute(h_iter_id);
    } while (h_iter_id != halfedge_id);
    m_fh_accessor.index_access().scalar_attribute(face_id) = halfedge_id;
}

void AtomicOperation::set_next(long halfedge_id, long next_halfedge_id)
{
    m_hn_accessor.index_access().scalar_attribute(halfedge_id) = next_halfedge_id;
    m_hp_accessor.index_access().scalar_attribute(next_halfedge_id) = halfedge_id;
}

void AtomicOperation::set_vertex(long halfedge_id, long vertex_id)
{
    long h_iter_id = halfedge_id;
    do {
        m_hv_accessor.index_access().scalar_attribute(h_iter_id) = vertex_id;
        h_iter_id = m_hp_accessor.index_access().scalar_attribute(implicit_opp(h_iter_id));
    } while (h_iter_id != halfedge_id);
    m_vh_accessor.index_access().scalar_attribute(vertex_id) = halfedge_id;
}

long AtomicOperation::get_face(long halfedge_id)
{
    return m_hf_accessor.index_access().scalar_attribute(halfedge_id);
}

long AtomicOperation::get_next(long halfedge_id)
{
    return m_hn_accessor.index_access().scalar_attribute(halfedge_id);
}

long AtomicOperation::get_vertex(long halfedge_id)
{
    return m_hv_accessor.index_access().scalar_attribute(halfedge_id);
}

long AtomicOperation::get_halfedge_from_tuple(const Tuple& t)
{
    return mesh().id(t, PrimitiveType::HalfEdge);
}

} // namespace wmtk::operations::polygon_mesh
