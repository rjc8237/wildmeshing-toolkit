#include "MakeEdge.hpp"

namespace wmtk::operations::polygon_mesh {

MakeEdge::MakeEdge(Mesh& m, const OperationSettings<MakeEdge>& settings)
    : AtomicOperation(m)
    , m_settings(settings)
{}

std::string MakeEdge::name() const
{
    return "polygon_mesh_make_edge";
}

bool MakeEdge::execute()
{
    long e_id = new_edge();
    long hij_id = 2 * e_id;
    long hji_id = 2 * e_id + 1;

    // Set new halfedges to be each others next
    set_next(hij_id, hji_id);
    set_next(hji_id, hij_id);

    // Set halfedges to have two distinct tips
    long vi_id = new_vertex();
    long vj_id = new_vertex();
    set_vertex(hij_id, vj_id);
    set_vertex(hji_id, vi_id);

    // Set halfedges to have a common (not hole) face
    long f_id = new_face(false);
    set_face(hij_id, f_id);
    set_face(hji_id, f_id);

    // Record output tuple
    m_output_tuple = mesh().tuple_from_id(PrimitiveType::Edge, e_id);

    return true;
}

Tuple MakeEdge::return_tuple() const
{
    return m_output_tuple;
}

} // namespace wmtk::operations::polygon_mesh
