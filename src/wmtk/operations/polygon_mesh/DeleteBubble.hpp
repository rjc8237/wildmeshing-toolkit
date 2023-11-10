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
};

namespace polygon_mesh {

/**
 * @class Atomic operation for removing a bubble component consisting of a single edge with two
 * distinct endpoints and a single face.
 */
class DeleteBubble : public AtomicOperation
{
public:
    // constructor for default factory pattern construction
    DeleteBubble(Mesh& m, const Tuple& e, const OperationSettings<DeleteBubble>& settings);

    std::string name() const override;

    /**
     * @brief Check the precondition that the tuple specifies a bubble component edge.
     *
     * @return true if the tuple specifies a bubble edge
     * @return false otherwise
     */
    bool precondition() const;

protected:
    bool execute() override;

private:
    const OperationSettings<DeleteBubble>& m_settings;
    Tuple m_tuple;
};

} // namespace polygon_mesh
} // namespace wmtk::operations
