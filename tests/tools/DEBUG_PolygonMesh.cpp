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
} // namespace wmtk::tests
