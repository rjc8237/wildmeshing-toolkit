#include "SplitVertex.hpp"
#include <wmtk/Simplex.hpp>
#include "FillHole.hpp"
#include "MakeEdge.hpp"
#include "MakeHole.hpp"
#include "Splice.hpp"

namespace wmtk::operations::polygon_mesh {

SplitVertex::SplitVertex(
    Mesh& m,
    const Tuple& h,
    const Tuple& g,
    const OperationSettings<SplitVertex>& settings)
    : PolygonMeshOperation(m)
    , m_first_tuple(h)
    , m_second_tuple(g)
    , m_settings(settings)
{}

std::string SplitVertex::name() const
{
    return "polygon_mesh_split_face";
}

bool SplitVertex::execute()
{
    assert(precondition());

    // Check if the faces adjacent to the halfedge are holes
    std::array<bool, 2> is_f_hole = {
        mesh().is_hole_face(m_first_tuple),
        mesh().is_hole_face(m_second_tuple)};

    // Splice the two halfedges together
    OperationSettings<Splice> splice_settings;
    if (!Splice(mesh(), m_first_tuple, m_second_tuple, splice_settings)()) {
        return false;
    }

    // Generate a new edge and get Tuples for the new halfedges
    OperationSettings<MakeEdge> make_edge_settings;
    MakeEdge make_edge(mesh(), make_edge_settings);
    if (!make_edge()) {
        return false;
    }
    Tuple e_tuple = make_edge.return_tuple();
    std::array<Tuple, 2> e_tuples = {e_tuple, mesh().opp_halfedge(e_tuple)};

    // The split is now done by simply splicing the new edge into the face
    if (!Splice(mesh(), e_tuples[0], m_second_tuple, splice_settings)()) {
        return false;
    }
    if (!Splice(mesh(), e_tuples[1], m_first_tuple, splice_settings)()) {
        return false;
    }

    // Mark faces as holes or filled as necessary
    for (long i = 0; i < 2; ++i) {
        if ((is_f_hole[i]) && (!mesh().is_hole_face(e_tuples[i]))) {
            OperationSettings<MakeHole> make_hole_settings;
            if (!MakeHole(mesh(), e_tuples[i], make_hole_settings)()) {
                return false;
            }
        } else if ((!is_f_hole[i]) && (mesh().is_hole_face(e_tuples[i]))) {
            OperationSettings<FillHole> fill_hole_settings;
            if (!FillHole(mesh(), e_tuples[i], fill_hole_settings)()) {
                return false;
            }
        }
    }

    // Set output tuple to first e tuple
    m_output_tuple = e_tuples[0];

    return true;
}

Tuple SplitVertex::return_tuple() const
{
    return m_output_tuple;
}

bool SplitVertex::precondition()
{
    // Halfedge tip vertices must be the same (switch vertex from default base to tip)
    Simplex v0(PrimitiveType::Vertex, mesh().switch_vertex(m_first_tuple));
    Simplex v1(PrimitiveType::Vertex, mesh().switch_vertex(m_second_tuple));
    if (!mesh().simplices_are_equal(v0, v1)) {
        return false;
    }

    // Halfedges cannot be the same
    Simplex h0(PrimitiveType::HalfEdge, m_first_tuple);
    Simplex h1(PrimitiveType::HalfEdge, m_second_tuple);
    if (mesh().simplices_are_equal(h0, h1)) {
        return false;
    }

    return true;
}

} // namespace wmtk::operations::polygon_mesh
