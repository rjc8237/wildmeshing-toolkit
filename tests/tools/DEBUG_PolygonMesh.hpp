#pragma once
#include <array>
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

    /**
     * @brief Find a ccw tuple in a given face originating from a given vertex.
     *
     * @param vid: id of a vertex in a face
     * @param fid: id of a face in the mesh
     * @return tuple of a ccw halfedge with a given face and base vertex
     */
    Tuple halfedge_tuple_from_vertex_in_face(long vid, long fid) const;

    /**
     * @brief Count the number of primitives in the mesh
     *
     * @return counts of the vertices, edges, faces, and halfedges respectively
     */
    std::array<long, 4> primitive_counts() const;

    /**
     * @brief Find per halfedge ids and counts of the connected components of the mesh.
     *
     * Hole faces are ignored and treated as interior faces in the computation.
     *
     * @param component_ids: per halfedge component ids
     * @return number of connected components
     */
    long find_simply_connected_components(std::vector<long>& component_ids) const;

    /**
     * @brief Count the number of connected components of the mesh.
     *
     * Hole faces are ignored and treated as interior faces in the computation.
     *
     * @return number of connected components
     */
    long count_simply_connected_components() const;

    /**
     * @brief Count the number of boundary loops in the mesh.
     *
     * @return number of boundary loops
     */
    long count_boundary_loops() const;

    /**
     * @brief Count the number of hole faces in the mesh
     *
     * @return number of hole faces
     */
    long count_hole_faces() const;

    /**
     * @brief Compute the Euler characteristic of the mesh
     *
     * @return Euler characteristic
     */
    long euler_characteristic() const;

    /**
     * @brief Compute the genus of a connected mesh
     *
     * @return genus of the mesh
     */
    long genus() const;
};

} // namespace wmtk::tests
