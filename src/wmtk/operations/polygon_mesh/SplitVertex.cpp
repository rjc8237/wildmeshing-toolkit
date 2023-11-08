#include "SplitVertex.hpp"
#include <wmtk/Simplex.hpp>
#include "MakeEdge.hpp"
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
    Tuple h0 = make_edge.return_tuple();
    Tuple h1 = mesh().opp_halfedge(h0);

    // The split is now done by simply splicing the new edge into the face
    if (!Splice(mesh(), h0, m_first_tuple, splice_settings)()) {
        return false;
    }
    if (!Splice(mesh(), h1, m_second_tuple, splice_settings)()) {
        return false;
    }

    return true;
}

bool SplitVertex::precondition()
{
    // Vertices must be the same
    Simplex f0(PrimitiveType::Vertex, m_first_tuple);
    Simplex f1(PrimitiveType::Vertex, m_second_tuple);
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
