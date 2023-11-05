#include "DeleteBubble.hpp"

namespace {
// Implicit opposite map for halfedges paired as [e] = {2*e, 2*e + 1}
long implicit_opp(long h)
{
    return ((h % 2) == 0) ? (h + 1) : (h - 1);
}
} // namespace

namespace wmtk::operations::polygon_mesh {
DeleteBubble::DeleteBubble(Mesh& m, const Tuple& e, const OperationSettings<DeleteBubble>& settings)
    : AtomicOperation(m)
    , m_settings(settings)
    , m_tuple(e)
{}

std::string DeleteBubble::name() const
{
    return "polygon_mesh_delete_bubble";
}

bool DeleteBubble::execute()
{
    assert(precondition());

    // Get element ids
    long hij_id = get_halfedge_from_tuple(m_tuple);
    long hji_id = implicit_opp(hij_id);
    long eij_id = hij_id / 2;
    long f_id = get_face(hij_id);
    long vj_id = get_vertex(hij_id);
    long vi_id = get_vertex(hji_id);


    // Delete elements
    delete_face(f_id);
    delete_vertex(vi_id);
    delete_vertex(vj_id);
    delete_edge(eij_id);

    return true;
}

bool DeleteBubble::precondition()
{
    // Get element ids
    long hij_id = get_halfedge_from_tuple(m_tuple);
    long hji_id = implicit_opp(hij_id);
    long f_id = get_face(hij_id);

    // Check the edge is actually an isolated bubble
    if ((get_next(hij_id) != hji_id) || (get_next(hji_id) != hij_id)) {
        return false;
    }
    if (get_face(hji_id) != f_id) {
        return false;
    }

    return true;
}

} // namespace wmtk::operations::polygon_mesh
