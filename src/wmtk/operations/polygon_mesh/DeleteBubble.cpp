#include "DeleteBubble.hpp"

namespace wmtk::operations::polygon_mesh {
DeleteBubble::DeleteBubble(Mesh& m, const Tuple& e, const OperationSettings<DeleteBubble>& settings)
    : AtomicOperation(m)
    , m_settings(settings)
{}

std::string DeleteBubble::name() const
{
    return "polygon_mesh_delete_bubble";
}

bool DeleteBubble::execute()
{
    // Get element ids
    long hij_id = get_halfedge_from_tuple(m_tuple);
    long hji_id = implicit_opp(hij_id);
    long eij_id = hij_id / 2;
    long f_id = get_face(hij_id);
    long vj_id = get_vertex(hij_id);
    long vi_id = get_vertex(hji_id);

    // Check the edge is actually an isolated bubble
    if ((get_next(hij_id) != hji_id) || (get_next(hji_id) != hij_id)) {
        return false;
    }
    if (get_face(hji_id) != f_id) {
        return false;
    }

    // Delete elements
    delete_face(f_id);
    delete_vertex(vi_id);
    delete_vertex(vj_id);
    delete_edge(eij_id);

    return true;
}

} // namespace wmtk::operations::polygon_mesh
