#include "PolygonMeshOperation.hpp"
#include <wmtk/PolygonMesh.hpp>

namespace wmtk::operations::polygon_mesh {
PolygonMeshOperation::PolygonMeshOperation(PolygonMesh& m)
    : m_mesh(m)
    , m_hash_accessor(get_hash_accessor(m))
{}

PolygonMeshOperation::PolygonMeshOperation(Mesh& m)
    : PolygonMeshOperation(dynamic_cast<PolygonMesh&>(m))
{}

Mesh& PolygonMeshOperation::base_mesh() const
{
    return m_mesh;
}
Accessor<long>& PolygonMeshOperation::hash_accessor()
{
    return m_hash_accessor;
}

PolygonMesh& PolygonMeshOperation::mesh() const
{
    return m_mesh;
}
} // namespace wmtk::operations::polygon_mesh
