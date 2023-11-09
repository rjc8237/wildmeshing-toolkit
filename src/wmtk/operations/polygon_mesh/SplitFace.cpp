#include "SplitFace.hpp"
#include <wmtk/Simplex.hpp>
#include "MakeEdge.hpp"
#include "MakeHole.hpp"
#include "Splice.hpp"

namespace wmtk::operations::polygon_mesh {

SplitFace::SplitFace(
    Mesh& m,
    const Tuple& h,
    const Tuple& g,
    const OperationSettings<SplitFace>& settings)
    : PolygonMeshOperation(m)
    , m_first_tuple(h)
    , m_second_tuple(g)
    , m_settings(settings)
{}

std::string SplitFace::name() const
{
    return "polygon_mesh_split_face";
}

bool SplitFace::execute()
{
    assert(precondition());

    // Generate a new edge and get Tuples for the new halfedges
    OperationSettings<MakeEdge> make_edge_settings;
    MakeEdge make_edge(mesh(), make_edge_settings);
    if (!make_edge()) {
        return false;
    }
    Tuple h0 = make_edge.return_tuple();
    Tuple h1 = mesh().opp_halfedge(h0);

    // Make the new bubble face a hole if the face to split is a hole
    if (mesh().is_hole_face(m_first_tuple)) {
        OperationSettings<MakeHole> make_hole_settings;
        if (!MakeHole(mesh(), h0, make_hole_settings)()) {
            return false;
        }
    }

    // The split is done by simply splicing the new edge into the face
    OperationSettings<Splice> splice_settings;
    if (!Splice(mesh(), m_second_tuple, h0, splice_settings)()) {
        return false;
    }
    if (!Splice(mesh(), m_first_tuple, h1, splice_settings)()) {
        return false;
    }

    // Set output tuple to h0
    m_output_tuple = h0;

    return true;
}

Tuple SplitFace::return_tuple() const
{
    return m_output_tuple;
}

bool SplitFace::precondition()
{
    // Faces must be the same
    Simplex f0(PrimitiveType::Face, m_first_tuple);
    Simplex f1(PrimitiveType::Face, m_second_tuple);
    if (!mesh().simplices_are_equal(f0, f1)) {
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
