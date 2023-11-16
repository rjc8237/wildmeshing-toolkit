#include "Splice.hpp"
#include <algorithm>
#include "FillHole.hpp"

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
    std::array<long, 2> old_h_ids = {
        get_halfedge_from_tuple(m_first_tuple),
        get_halfedge_from_tuple(m_second_tuple)};
    std::array<long, 2> old_hn_ids = {get_next(old_h_ids[0]), get_next(old_h_ids[1])};

    // Do nothing for identical halfedges
    if (old_hn_ids[0] == old_hn_ids[1]) {
        return true;
    }

    // Old face and vertex indices
    std::array<long, 2> old_hf_ids = {get_face(old_h_ids[0]), get_face(old_h_ids[1])};
    std::array<long, 2> old_hv_ids = {get_vertex(old_h_ids[0]), get_vertex(old_h_ids[1])};

    // Swap next(h0) and next(h1) (with corresponding change in prev)
    set_next(old_h_ids[0], old_hn_ids[1]);
    set_next(old_h_ids[1], old_hn_ids[0]);

    // The new face is a hole iff both input faces are holes
    bool new_face_is_hole = (is_hole(old_hf_ids[0]) && is_hole(old_hf_ids[1]));

    // If the halfedges share a face, then it is split into two distinct faces
    if (old_hf_ids[0] == old_hf_ids[1]) {
        for (const auto& old_h_id : old_h_ids) {
            long new_f_id = new_face(new_face_is_hole);
            set_face(old_h_id, new_f_id);
        }
    }
    // If the halfedges have distinct faces, then they are merged into a single face
    else {
        long new_f_id = new_face(new_face_is_hole);
        for (const auto& old_h_id : old_h_ids) {
            set_face(old_h_id, new_f_id);
        }
    }
    // Delete old face ids
    for (const auto& old_hf_id : old_hf_ids) {
        delete_face(old_hf_id);
    }

    // If the halfedges share a tip vertex, then it is split into two distinct vertices
    if (old_hv_ids[0] == old_hv_ids[1]) {
        for (const auto& old_h_id : old_h_ids) {
            long new_v_id = new_vertex();
            set_vertex(old_h_id, new_v_id);
        }
    }
    // If the halfedges have distinct tips, then they are merged into a single vertex
    else {
        long new_v_id = new_vertex();
        for (const auto& old_h_id : old_h_ids) {
            set_vertex(old_h_id, new_v_id);
        }
    }
    // Delete old vertex ids
    for (const auto& old_hv_id : old_hv_ids) {
        delete_vertex(old_hv_id);
    }

    return true;
}

long Splice::next_after_splice(long h) const
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

bool Splice::is_hole_after_splice(long fid) const
{
    std::array<long, 2> old_h_ids = {
        get_halfedge_from_tuple(m_first_tuple),
        get_halfedge_from_tuple(m_second_tuple)};
    std::array<long, 2> old_hf_ids = {get_face(old_h_ids[0]), get_face(old_h_ids[1])};

    // Faces adjacent to spliced halfedges are holes iff both input are holes
    if ((fid == old_hf_ids[0]) || (fid == old_hf_ids[1])) {
        return (is_hole(old_hf_ids[0]) && is_hole(old_hf_ids[1]));
    } else {
        return is_hole(fid);
    }
}

bool Splice::precondition() const
{
    // Get old halfedge indices
    std::array<long, 2> old_h_ids = {
        get_halfedge_from_tuple(m_first_tuple),
        get_halfedge_from_tuple(m_second_tuple)};
    std::array<long, 2> old_hf_ids = {get_face(old_h_ids[0]), get_face(old_h_ids[1])};
    std::array<long, 2> old_hv_ids = {get_vertex(old_h_ids[0]), get_vertex(old_h_ids[1])};

    // Check boundary validity at each vertex of the face if merging a hole and interior face
    if (is_hole(old_hf_ids[0]) != is_hole(old_hf_ids[1])) {
        long h_iter = old_h_ids[0];
        do {
            if (!precondition_at_vertex(h_iter)) {
                return false;
            }
            h_iter = next_after_splice(h_iter);
        } while (h_iter != old_h_ids[0]);
    }
    // Check boundary validity at just the new vertex if merging distinct vertices
    else if (old_hv_ids[0] != old_hv_ids[1]) {
        if (!precondition_at_vertex(old_h_ids[0])) {
            return false;
        }
    }

    return true;
}


bool Splice::precondition_at_vertex(long halfedge_id) const
{
    // Get old halfedge indices
    std::array<long, 2> old_h_ids = {
        get_halfedge_from_tuple(m_first_tuple),
        get_halfedge_from_tuple(m_second_tuple)};
    std::array<long, 2> old_hf_ids = {get_face(old_h_ids[0]), get_face(old_h_ids[1])};

    // Check boundary validity at merged vertex
    long h_iter = halfedge_id;
    long num_boundary_edges = 0;
    bool prev_is_hole = is_hole_after_splice(get_face(h_iter));
    // Check if the current face is a hole
    do {
        // Circulate around the vertex
        h_iter = implicit_opp(next_after_splice(h_iter));

        // Check if switched from hole to interior or reverse
        bool curr_is_hole = is_hole_after_splice(get_face(h_iter));
        if (curr_is_hole != prev_is_hole) {
            num_boundary_edges++;
        }

        // Update current hole status
        prev_is_hole = curr_is_hole;
    } while (h_iter != halfedge_id);

    assert((num_boundary_edges % 2) == 0);
    return (num_boundary_edges <= 2);
}

bool Splice::is_topology_preserving() const
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
