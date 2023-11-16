#include "JoinVertex.hpp"
#include <wmtk/Simplex.hpp>
#include "DeleteBubble.hpp"
#include "FillHole.hpp"
#include "MakeHole.hpp"
#include "Splice.hpp"

namespace wmtk::operations::polygon_mesh {

JoinVertex::JoinVertex(Mesh& m, const Tuple& h, const OperationSettings<JoinVertex>& settings)
    : PolygonMeshOperation(m)
    , m_tuple(h)
    , m_settings(settings)
{}

std::string JoinVertex::name() const
{
    return "polygon_mesh_join_vertex";
}

bool JoinVertex::execute()
{
    assert(precondition());

    // Get relevant halfedge tuples
    std::array<Tuple, 2> h_tuples = {m_tuple, mesh().opp_halfedge(m_tuple)};
    std::array<Tuple, 2> hp_tuples = {
        mesh().prev_halfedge(h_tuples[0]),
        mesh().prev_halfedge(h_tuples[1])};

    // Check if the faces adjacent to the halfedge are holes
    std::array<bool, 2> is_f_hole = {
        mesh().is_hole_face(h_tuples[0]),
        mesh().is_hole_face(h_tuples[1])};

    // The join is done by simply splicing out the edge
    OperationSettings<Splice> splice_settings;
    if (!Splice(mesh(), hp_tuples[1], h_tuples[0], splice_settings)()) {
        return false;
    }
    if (!Splice(mesh(), hp_tuples[0], h_tuples[1], splice_settings)()) {
        return false;
    }
    if (!Splice(mesh(), hp_tuples[0], hp_tuples[1], splice_settings)()) {
        return false;
    }

    // Delete the spliced out edge bubble
    OperationSettings<DeleteBubble> delete_bubble_settings;
    DeleteBubble delete_bubble(mesh(), h_tuples[0], delete_bubble_settings);
    if (!delete_bubble()) {
        return false;
    }

    // Mark faces as holes or filled as necessary
    for (long i = 0; i < 2; ++i) {
        if ((is_f_hole[i]) && (!mesh().is_hole_face(hp_tuples[i]))) {
            OperationSettings<MakeHole> make_hole_settings;
            if (!MakeHole(mesh(), hp_tuples[i], make_hole_settings)()) {
                return false;
            }
        } else if ((!is_f_hole[i]) && (mesh().is_hole_face(hp_tuples[i]))) {
            OperationSettings<FillHole> fill_hole_settings;
            if (!FillHole(mesh(), hp_tuples[i], fill_hole_settings)()) {
                return false;
            }
        }
    }

    // Set the output tuple to the next halfedge after g0
    m_output_tuple = mesh().next_halfedge(hp_tuples[0]);

    return true;
}

Tuple JoinVertex::return_tuple() const
{
    return m_output_tuple;
}

bool JoinVertex::precondition()
{
    // If the edge is not a boundary, both vertices cannot be boundary
    Tuple v0_tuple = m_tuple;
    Tuple v1_tuple = mesh().opp_halfedge(m_tuple);
    bool v0_is_boundary = (mesh().is_boundary_vertex(v0_tuple));
    bool v1_is_boundary = (mesh().is_boundary_vertex(v1_tuple));
    bool e_is_boundary = (mesh().is_boundary_edge(m_tuple));
    if ((v0_is_boundary) && (v1_is_boundary) && (!e_is_boundary)) {
        return false;
    }

    // The previous halfedges of the edge must be distinct
    Simplex h0(PrimitiveType::HalfEdge, m_tuple);
    Simplex h1(PrimitiveType::HalfEdge, mesh().opp_halfedge(m_tuple));
    Simplex g0(PrimitiveType::HalfEdge, mesh().prev_halfedge(m_tuple));
    Simplex g1(PrimitiveType::HalfEdge, mesh().prev_halfedge(mesh().opp_halfedge(m_tuple)));
    if ((mesh().simplices_are_equal(h0, g0)) || (mesh().simplices_are_equal(h0, g1)) ||
        (mesh().simplices_are_equal(h1, g0)) || (mesh().simplices_are_equal(h1, g1))) {
        return false;
    }

    return true;
}

} // namespace wmtk::operations::polygon_mesh
