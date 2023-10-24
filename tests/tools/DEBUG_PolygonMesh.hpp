#pragma once
#include <wmtk/PolygonMesh.hpp>

namespace wmtk::tests {
class DEBUG_PolygonMesh : public PolygonMesh //, public virtual DEBUG_Mesh
{
public:
    using PolygonMesh::PolygonMesh;
    DEBUG_PolygonMesh(const PolygonMesh& m);
    DEBUG_PolygonMesh(PolygonMesh&& m);
    using PolygonMesh::operator=;


    bool operator==(const DEBUG_PolygonMesh& o) const;
    bool operator!=(const DEBUG_PolygonMesh& o) const;

    long id(const Tuple& tuple, PrimitiveType type) const override;
    using PolygonMesh::tuple_from_id;
};

} // namespace wmtk::tests
