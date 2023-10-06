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

    /**
     * @brief Generate a new vertex handle
     * 
     * TODO: Check if should be implemented using request_simplex_indices
     */
    long new_vertex();

    /**
     * @brief Generate a clone of an existing vertex handle
     * 
     * TODO: Check if this is necessary if handles are just indices
     */
    long clone_vertex(long v);

    /**
     * @brief Delete the vertex with the given id
     */
    void delete_vertex(long v);

    /**
     * @brief Set the id of the vertex corresponding to the given tuple to v
     */
    void set_vertex(const Tuple& tuple, long v);

    /**
     * @brief Generate a new face handle
     */
    long new_face();

    /**
     * @brief Generate a clone of an existing face handle
     */
    long clone_face(long f);

    /**
     * @brief Delete the face with the given id
     */
    void delete_face(long f);

    /**
     * @brief Set the id of the face corresponding to the given tuple to f
     */
    void set_face(const Tuple& tuple, long f);

    /**
     * @brief Generate a new edge handle
     */
    long new_edge();
private:
    PolygonMesh& m_mesh;
    Accessor<long> m_hash_accessor;
};


} // namespace operations::polygon_mesh
} // namespace wmtk
