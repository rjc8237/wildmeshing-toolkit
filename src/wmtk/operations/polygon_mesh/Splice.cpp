#include "Splice.hpp"


namespace {
// Implicit opposite map for halfedges paired as [e] = {2*e, 2*e + 1}
long implicit_opp(long h)
{
    return ((h % 2) == 0) ? (h + 1) : (h - 1);
}
} // namespace

namespace wmtk::operations::polygon_mesh {

Splice::Splice(Mesh& m, const Tuple& h, const Tuple& g, const OperationSettings<Splice>& settings)
    : AtomicOperation(m)
    , m_first_tuple(h)
    , m_second_tuple(g)
    , m_settings(settings)
{}

std::string Splice::name() const
{
    return "polygon_mesh_splice";
}

bool Splice::execute()
{
    assert(precondition());

    // Old halfedge indices
    long old_h_id = get_halfedge_from_tuple(m_first_tuple);
    long old_g_id = get_halfedge_from_tuple(m_second_tuple);
    long old_hn_id = get_next(old_h_id);
    long old_gn_id = get_next(old_g_id);

    // Do nothing for identical halfedges
    if (old_hn_id == old_gn_id) {
        return true;
    }

    // Old face indices
    long old_hf_id = get_face(old_h_id);
    long old_gf_id = get_face(old_g_id);

    // Old vertex indices
    long old_hv_id = get_vertex(old_h_id);
    long old_gv_id = get_vertex(old_g_id);

    // Swap next(h) and next(g) (with corresponding change in prev)
    set_next(old_h_id, old_gn_id);
    set_next(old_g_id, old_hn_id);

    // If the halfedges share a face, then it is split
    if (old_hf_id == old_gf_id) {
        bool split_face_is_hole = (is_hole(old_hf_id) && is_hole(old_gf_id));
        long new_hf_id = new_face(split_face_is_hole);
        long new_gf_id = new_face(split_face_is_hole);
        set_face(old_h_id, new_hf_id);
        set_face(old_g_id, new_gf_id);
    }
    // If the halfedges have distinct faces, then they are merged
    else {
        // If one of the faces is a hole, close it in the splice
        bool merged_face_is_hole = (is_hole(old_hf_id) && is_hole(old_gf_id));
        long new_f_id = new_face(merged_face_is_hole);
        set_face(old_h_id, new_f_id);
        set_face(old_g_id, new_f_id);
    }
    delete_face(old_hf_id);
    delete_face(old_gf_id);

    // If the halfedges share a vertex, then it is split
    if (old_hv_id == old_gv_id) {
        long new_hv_id = new_vertex();
        long new_gv_id = new_vertex();
        set_vertex(old_h_id, new_hv_id);
        set_vertex(old_g_id, new_gv_id);
    }
    // If the halfedges have distinct vertices, then they are merged
    else {
        long new_v_id = new_vertex();
        set_vertex(old_h_id, new_v_id);
        set_vertex(old_g_id, new_v_id);
    }
    delete_vertex(old_hv_id);
    delete_vertex(old_gv_id);

    return true;
}

long Splice::next_after_splice(long h)
{
    long old_h_id = get_halfedge_from_tuple(m_first_tuple);
    long old_g_id = get_halfedge_from_tuple(m_second_tuple);
    if (h == old_h_id) {
        return get_next(old_g_id);
    } else if (h == old_g_id) {
        return get_next(old_h_id);
    } else {
        return get_next(h);
    }
}

bool Splice::precondition()
{
    long old_h_id = get_halfedge_from_tuple(m_first_tuple);
    long old_g_id = get_halfedge_from_tuple(m_second_tuple);
    long old_hv_id = get_vertex(old_h_id);
    long old_gv_id = get_vertex(old_g_id);

    // Splitting a vertex is always manifold
    if (old_hv_id == old_gv_id) {
        return true;
    }

    // Check boundary validity
    long h_iter = old_h_id;
    long num_boundary_edges = 0;
    bool prev_is_hole = is_hole(get_face(old_h_id));
    do {
        // Circulate around the vertex
        h_iter = implicit_opp(next_after_splice(h_iter));

        // Check if switched from hole to interior or reverse
        long curr_is_hole = is_hole(get_face(h_iter));
        if (curr_is_hole != prev_is_hole) {
            num_boundary_edges++;
        }

        // Update current hole status
        prev_is_hole = curr_is_hole;
    } while (h_iter != old_h_id);

    assert((num_boundary_edges % 2) == 0);
    return (num_boundary_edges <= 2);
}

bool Splice::is_topology_preserving()
{
    // halfedge indices
    long h_id = get_halfedge_from_tuple(m_first_tuple);
    long g_id = get_halfedge_from_tuple(m_second_tuple);

    // face indices
    long hf_id = get_face(h_id);
    long gf_id = get_face(g_id);

    // Old vertex indices
    long hv_id = get_vertex(h_id);
    long gv_id = get_vertex(g_id);

    return ((hf_id == gf_id) && (hv_id != gv_id)) || ((hf_id != gf_id) && (hv_id == gv_id));
}


} // namespace wmtk::operations::polygon_mesh
