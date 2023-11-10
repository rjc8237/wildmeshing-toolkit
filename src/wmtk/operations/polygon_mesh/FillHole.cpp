#include "FillHole.hpp"

namespace {
// Implicit opposite map for halfedges paired as [e] = {2*e, 2*e + 1}
long implicit_opp(long h)
{
    return ((h % 2) == 0) ? (h + 1) : (h - 1);
}
} // namespace

namespace wmtk::operations::polygon_mesh {

FillHole::FillHole(Mesh& m, const Tuple& f, const OperationSettings<FillHole>& settings)
    : AtomicOperation(m)
    , m_settings(settings)
    , m_tuple(f)
{}

std::string FillHole::name() const
{
    return "polygon_mesh_fill_hole";
}

bool FillHole::execute()
{
    assert(precondition());

    fill_hole(get_face(get_halfedge_from_tuple(m_tuple)));

    return true;
}

bool FillHole::precondition() const
{
    long halfedge_id = get_halfedge_from_tuple(m_tuple);

    // Nothing to check if the face is already full
    if (!is_hole(get_face(halfedge_id))) {
        return true;
    }

    // Check boundary validity at each vertex of the face
    long h_iter = halfedge_id;
    do {
        if (!precondition_at_vertex(h_iter)) {
            return false;
        }
        h_iter = get_next(h_iter);
    } while (h_iter != halfedge_id);

    return true;
}

bool FillHole::precondition_at_vertex(long halfedge_id) const
{
    // Check boundary validity
    long h_iter = halfedge_id;
    long num_boundary_edges = 0;
    long new_hole_face = get_face(get_halfedge_from_tuple(m_tuple));
    bool prev_is_hole = (is_hole(get_face(h_iter)) && (get_face(h_iter) != new_hole_face));
    do {
        // Circulate around the vertex
        h_iter = implicit_opp(get_next(h_iter));

        // Check if switched from hole to interior or reverse
        long curr_is_hole = (is_hole(get_face(h_iter)) && (get_face(h_iter) != new_hole_face));
        if (curr_is_hole != prev_is_hole) {
            num_boundary_edges++;
        }

        // Update current hole status
        prev_is_hole = curr_is_hole;
    } while (h_iter != halfedge_id);

    assert((num_boundary_edges % 2) == 0);
    return (num_boundary_edges <= 2);
}

} // namespace wmtk::operations::polygon_mesh
