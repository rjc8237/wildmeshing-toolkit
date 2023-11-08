#include "JoinVertex.hpp"
#include <wmtk/Simplex.hpp>
#include "DeleteBubble.hpp"
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

    // The join is done by simply splicing out the edge
    Tuple h0_tuple = m_tuple;
    Tuple h1_tuple = mesh().opp_halfedge(h0_tuple);
    Tuple g0_tuple = mesh().prev_halfedge(h0_tuple);
    Tuple g1_tuple = mesh().prev_halfedge(h1_tuple);
    OperationSettings<Splice> splice_settings;
    if (!Splice(mesh(), g1_tuple, h0_tuple, splice_settings)()) {
        return false;
    }
    if (!Splice(mesh(), g0_tuple, h1_tuple, splice_settings)()) {
        return false;
    }
    if (!Splice(mesh(), g0_tuple, g1_tuple, splice_settings)()) {
        return false;
    }

    // Delete the spliced out edge bubble
    OperationSettings<DeleteBubble> delete_bubble_settings;
    DeleteBubble delete_bubble(mesh(), h0_tuple, delete_bubble_settings);
    if (!delete_bubble()) {
        return false;
    }

    return true;
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

    // At least one vertex must have valence greater than 1
    Simplex h0(PrimitiveType::HalfEdge, m_tuple);
    Simplex h1(PrimitiveType::HalfEdge, mesh().opp_halfedge(m_tuple));
    Simplex g0(PrimitiveType::HalfEdge, mesh().next_halfedge(m_tuple));
    Simplex g1(PrimitiveType::HalfEdge, mesh().next_halfedge(mesh().opp_halfedge(m_tuple)));
    if ((mesh().simplices_are_equal(g0, h1)) && (mesh().simplices_are_equal(g1, h0))) {
        return false;
    }

    return true;
}

} // namespace wmtk::operations::polygon_mesh
