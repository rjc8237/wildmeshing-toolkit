#include "JoinFace.hpp"
#include <wmtk/Simplex.hpp>
#include "DeleteBubble.hpp"
#include "MakeHole.hpp"
#include "Splice.hpp"

namespace wmtk::operations::polygon_mesh {

JoinFace::JoinFace(Mesh& m, const Tuple& h, const OperationSettings<JoinFace>& settings)
    : PolygonMeshOperation(m)
    , m_tuple(h)
    , m_settings(settings)
{}

std::string JoinFace::name() const
{
    return "polygon_mesh_join_face";
}

bool JoinFace::execute()
{
    assert(precondition());

    // Determine if the new face is a hole
    Tuple h0_tuple = m_tuple;
    Tuple h1_tuple = mesh().opp_halfedge(h0_tuple);
    bool is_hole = ((mesh().is_hole_face(h0_tuple)) && (mesh().is_hole_face(h1_tuple)));

    // The join is done by simply splicing out the edge
    OperationSettings<Splice> splice_settings;
    Tuple g = mesh().prev_halfedge(h0_tuple);
    if (!Splice(mesh(), mesh().prev_halfedge(h1_tuple), h0_tuple, splice_settings)()) {
        return false;
    }
    if (!Splice(mesh(), mesh().prev_halfedge(h0_tuple), h1_tuple, splice_settings)()) {
        return false;
    }

    // Delete the spliced out edge bubble
    OperationSettings<DeleteBubble> delete_bubble_settings;
    DeleteBubble delete_bubble(mesh(), h0_tuple, delete_bubble_settings);
    if (!delete_bubble()) {
        return false;
    }

    // Mark the new face as a hole
    if (is_hole) {
        OperationSettings<MakeHole> make_hole_settings;
        MakeHole make_hole(mesh(), g, make_hole_settings);
    }

    return true;
}

bool JoinFace::precondition()
{
    // Faces of the edge must be distinct
    Simplex f0(PrimitiveType::Face, m_tuple);
    Simplex f1(PrimitiveType::Face, mesh().opp_halfedge(m_tuple));
    return (!mesh().simplices_are_equal(f0, f1));
}

} // namespace wmtk::operations::polygon_mesh
