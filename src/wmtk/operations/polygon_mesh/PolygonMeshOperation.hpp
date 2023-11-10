#pragma once
#include <wmtk/operations/Operation.hpp>

namespace wmtk {
class PolygonMesh;

namespace operations::polygon_mesh {

class PolygonMeshOperation : public virtual Operation
{
public:
    PolygonMeshOperation(PolygonMesh& m);
    // internally will try dynamic casting to check for mistakes
    PolygonMeshOperation(Mesh& m);

protected:
    PolygonMesh& mesh() const;
    Mesh& base_mesh() const override;
    Accessor<long>& hash_accessor() override;

private:
    PolygonMesh& m_mesh;
    Accessor<long> m_hash_accessor;
};


} // namespace operations::polygon_mesh
} // namespace wmtk
