#include "DEBUG_PolygonMesh.hpp"

namespace wmtk::tests {

DEBUG_PolygonMesh::DEBUG_PolygonMesh(const PolygonMesh& m)
    : PolygonMesh(m)
{}
DEBUG_PolygonMesh::DEBUG_PolygonMesh(PolygonMesh&& m)
    : PolygonMesh(std::move(m))
{}


bool DEBUG_PolygonMesh::operator==(const DEBUG_PolygonMesh& o) const
{
    return static_cast<const PolygonMesh&>(*this) == static_cast<const PolygonMesh&>(o);
}
bool DEBUG_PolygonMesh::operator!=(const DEBUG_PolygonMesh& o) const
{
    return !(*this == o);
}

long DEBUG_PolygonMesh::id(const Tuple& tuple, PrimitiveType type) const
{
    return PolygonMesh::id(tuple, type);
}

Tuple DEBUG_PolygonMesh::halfedge_tuple_from_vertex_in_face(long vid, long fid) const
{
    Tuple f_tuple = tuple_from_id(PrimitiveType::Face, fid);
    Tuple f_tuple_iter = f_tuple;
    do {
        if (id(f_tuple_iter, PrimitiveType::Vertex) == vid) {
            break;
        }
        f_tuple_iter = next_halfedge(f_tuple_iter);
    } while (id(f_tuple_iter, PrimitiveType::HalfEdge) != id(f_tuple, PrimitiveType::HalfEdge));

    assert(id(f_tuple_iter, PrimitiveType::Vertex) == vid);
    return f_tuple_iter;
}
} // namespace wmtk::tests
