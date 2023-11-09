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
    , m_next_accessor(m.create_accessor<long>(m.m_next_handle))
    , m_prev_accessor(m.create_accessor<long>(m.m_prev_handle))
    , m_to_accessor(m.create_accessor<long>(m.m_to_handle))
    , m_out_accessor(m.create_accessor<long>(m.m_out_handle))
    , m_hf_accessor(m.create_accessor<long>(m.m_hf_handle))
    , m_fh_accessor(m.create_accessor<long>(m.m_fh_handle))
    , m_f_is_hole_accessor(m.create_accessor<char>(m.m_f_is_hole_handle))
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
    long fid = mesh().request_simplex_indices(PrimitiveType::Face, 1)[0];
    mesh().get_index_access(m_f_is_hole_accessor).scalar_attribute(fid) = is_hole ? 1 : 0;
    return fid;
}

long AtomicOperation::new_edge()
{
    // Get new edge id
    mesh().reserve_attributes(PrimitiveType::Edge, mesh().capacity(PrimitiveType::Edge) + 1);
    long eid = mesh().request_simplex_indices(PrimitiveType::Edge, 1)[0];

    // Also make attributes for corresponding halfedge ids
    size_t primitive_id = get_primitive_type_id(PrimitiveType::HalfEdge);
    long new_halfedge_capacity = mesh().capacity(PrimitiveType::HalfEdge) + 2;
    mesh().reserve_attributes(PrimitiveType::HalfEdge, new_halfedge_capacity);
    mesh().m_attribute_manager.m_capacities[primitive_id] = new_halfedge_capacity;
    mesh().get_index_access(m_h_flag_accessor).scalar_attribute(2 * eid) |= 0x1;
    mesh().get_index_access(m_h_flag_accessor).scalar_attribute(2 * eid + 1) |= 0x1;

    return eid;
}

long AtomicOperation::new_vertex()
{
    mesh().reserve_attributes(PrimitiveType::Vertex, mesh().capacity(PrimitiveType::Vertex) + 1);
    return mesh().request_simplex_indices(PrimitiveType::Vertex, 1)[0];
}

void AtomicOperation::delete_face(long face_id)
{
    mesh().get_index_access(m_f_flag_accessor).scalar_attribute(face_id) = 0;
}

void AtomicOperation::delete_edge(long edge_id)
{
    mesh().get_index_access(m_e_flag_accessor).scalar_attribute(edge_id) = 0;
    mesh().get_index_access(m_h_flag_accessor).scalar_attribute(2 * edge_id) = 0;
    mesh().get_index_access(m_h_flag_accessor).scalar_attribute(2 * edge_id + 1) = 0;
}

void AtomicOperation::delete_vertex(long vertex_id)
{
    mesh().get_index_access(m_v_flag_accessor).scalar_attribute(vertex_id) = 0;
}

void AtomicOperation::set_face(long halfedge_id, long face_id)
{
    long h_iter_id = halfedge_id;
    do {
        mesh().get_index_access(m_hf_accessor).scalar_attribute(h_iter_id) = face_id;
        h_iter_id = mesh().get_index_access(m_next_accessor).scalar_attribute(h_iter_id);
    } while (h_iter_id != halfedge_id);
    mesh().get_index_access(m_fh_accessor).scalar_attribute(face_id) = halfedge_id;
}

void AtomicOperation::set_next(long halfedge_id, long next_halfedge_id)
{
    mesh().get_index_access(m_next_accessor).scalar_attribute(halfedge_id) = next_halfedge_id;
    mesh().get_index_access(m_prev_accessor).scalar_attribute(next_halfedge_id) = halfedge_id;
}

void AtomicOperation::set_vertex(long halfedge_id, long vertex_id)
{
    long h_iter_id = halfedge_id;
    do {
        mesh().get_index_access(m_to_accessor).scalar_attribute(h_iter_id) = vertex_id;
        h_iter_id =
            mesh().get_index_access(m_prev_accessor).scalar_attribute(implicit_opp(h_iter_id));
    } while (h_iter_id != halfedge_id);
    mesh().get_index_access(m_out_accessor).scalar_attribute(vertex_id) = implicit_opp(halfedge_id);
}


bool AtomicOperation::is_hole(long face_id) const
{
    return mesh().is_hole_face(mesh().tuple_from_id(PrimitiveType::Face, face_id));
}

void AtomicOperation::make_hole(long face_id)
{
    mesh().get_index_access(m_f_is_hole_accessor).scalar_attribute(face_id) |= 0x1;
}

void AtomicOperation::fill_hole(long face_id)
{
    mesh().get_index_access(m_f_is_hole_accessor).scalar_attribute(face_id) = 0;
}

long AtomicOperation::get_face(long halfedge_id) const
{
    return mesh().get_index_access(m_hf_accessor).scalar_attribute(halfedge_id);
}

long AtomicOperation::get_next(long halfedge_id) const
{
    return mesh().get_index_access(m_next_accessor).scalar_attribute(halfedge_id);
}

long AtomicOperation::get_vertex(long halfedge_id) const
{
    return mesh().get_index_access(m_to_accessor).scalar_attribute(halfedge_id);
}

long AtomicOperation::get_halfedge_from_tuple(const Tuple& t) const
{
    return mesh().id(t, PrimitiveType::HalfEdge);
}

} // namespace wmtk::operations::polygon_mesh
