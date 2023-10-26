#pragma once
#include <wmtk/operations/TupleOperation.hpp>
#include "AtomicOperation.hpp"
#include "PolygonMeshOperation.hpp"

namespace wmtk::operations {
namespace polygon_mesh {
class DeleteBubble;
}

template <>
struct OperationSettings<polygon_mesh::DeleteBubble>
{
    // TODO Determine if operation settings are necessary
};

namespace polygon_mesh {

/**
 * @class Atomic operation for removing a bubble component with a single edge.
 */
class DeleteBubble : public AtomicOperation
{
public:
    // constructor for default factory pattern construction
    DeleteBubble(Mesh& m, const Tuple& e, const OperationSettings<DeleteBubble>& settings);

    std::string name() const override;

    using PolygonMeshOperation::hash_accessor;

protected:
    bool execute() override;

private:
    const OperationSettings<DeleteBubble>& m_settings;
};

} // namespace polygon_mesh
} // namespace wmtk::operations
