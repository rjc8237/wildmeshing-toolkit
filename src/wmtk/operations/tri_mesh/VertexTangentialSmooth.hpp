#pragma once
#include <optional>
#include <wmtk/invariants/InvariantCollection.hpp>
#include <wmtk/operations/TupleOperation.hpp>
#include "TriMeshOperation.hpp"
#include "VertexAttributesUpdateBase.hpp"
#include "VertexLaplacianSmooth.hpp"

namespace wmtk::operations {
namespace tri_mesh {
class VertexTangentialSmooth;
}

template <>
struct OperationSettings<tri_mesh::VertexTangentialSmooth>
{
    OperationSettings<tri_mesh::VertexLaplacianSmooth> smooth_settings;

    double damping_factor = 1.0;
};

namespace tri_mesh {
class VertexTangentialSmooth : public VertexLaplacianSmooth
{
public:
    VertexTangentialSmooth(
        Mesh& m,
        const Tuple& t,
        const OperationSettings<VertexTangentialSmooth>& settings);

    std::string name() const override;

    static PrimitiveType primitive_type() { return PrimitiveType::Vertex; }

protected:
    bool before() const override;
    bool execute() override;

private:
    const OperationSettings<VertexTangentialSmooth>& m_settings;
};

} // namespace tri_mesh
} // namespace wmtk::operations
