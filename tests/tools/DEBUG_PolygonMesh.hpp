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

    template <typename T>
    attribute::AccessorBase<T> create_base_accessor(const MeshAttributeHandle<T>& handle)
    {
        return attribute::AccessorBase<T>(*this, handle);
    }

    template <typename T>
    attribute::AccessorBase<T> create_const_base_accessor(
        const MeshAttributeHandle<T>& handle) const
    {
        return attribute::AccessorBase<T>(const_cast<DEBUG_PolygonMesh&>(*this), handle);
    }
    template <typename T>
    attribute::AccessorBase<T> create_base_accessor(const MeshAttributeHandle<T>& handle) const
    {
        return create_const_base_accessor(handle);
    }

    const MeshAttributeHandle<long>& next_handle() const;

    /**
     * @brief Perform a complete validity check to ensure the polygon mesh data represents
     * a valid manifold surface.
     *
     * @return true if the connectivity is valid
     * @return false otherwise
     */
    bool is_connectivity_valid() const override;

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
     * @brief Determine if a vertex is a hole vertex (i.e., surrounded by hole faces)
     *
     * @param v: tuple representing a vertex
     * @return true if v is a hole vertex
     * @return false otherwise
     */
    bool is_hole_vertex(const Tuple& v) const;

    /**
     * @brief Determine if an edge is a hole edge (i.e., surrounded by hole faces)
     *
     * @param e: tuple representing a edge
     * @return true if e is a hole edge
     * @return false otherwise
     */
    bool is_hole_edge(const Tuple& e) const;

    /**
     * @brief Determine if a primitive specified by a tuple is a hole.
     *
     * @param tuple: tuple specifying a primitive
     * @param type: type of the primitive
     * @return true if the primitive is a hole
     * @return false otherwise
     */
    bool is_hole(const Tuple& tuple, PrimitiveType type) const;

    /**
     * @brief Count the number of primitives of a given type that are holes.
     *
     * @param type: type of hole primitives to count
     * @return long number of hole primitives
     */
    long count_hole_primitives(PrimitiveType type) const;

    /**
     * @brief Count the number of hole faces in the mesh
     *
     * @return number of hole faces
     */
    long count_hole_faces() const;

    /**
     * @brief Count the number of interior faces in the mesh
     *
     * @return number of interior faces
     */
    long count_interior_faces() const;

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

    /**
     * @brief Determine if the local connectivity for a primitive specified by tuple is valid
     *
     * @param tuple: tuple to check
     * @param type: primitive of the tuple to check
     * @return true iff the connectivity corresponding to the primitive is valid
     */
    bool is_local_connectivity_valid(const Tuple& tuple, PrimitiveType type) const;

private:
    /**
     * @brief Determine if the local connectivity for a vertex with a given id is valid
     *
     * @param vid: id of the vertex
     * @return true iff the id specifies a valid vertex and the connectivity involving it is valid
     */
    bool is_vertex_connectivity_valid(long vid) const;

    /**
     * @brief Determine if the local connectivity for a edge with a given id is valid
     *
     * @param eid: id of the edge
     * @return true iff the id specifies a valid edge and the connectivity involving it is valid
     */
    bool is_edge_connectivity_valid(long eid) const;

    /**
     * @brief Determine if the local connectivity for a face with a given id is valid
     *
     * @param fid: id of the face
     * @return true iff the id specifies a valid face and the connectivity involving it is valid
     */
    bool is_face_connectivity_valid(long fid) const;

    /**
     * @brief Determine if the local connectivity for a halfedge with a given id is valid
     *
     * @param hid: id of the halfedge
     * @return true iff the id specifies a valid halfedge and the connectivity involving it is valid
     */
    bool is_halfedge_connectivity_valid(long hid) const;

    /**
     * @brief Count the number of vertices in the mesh (including hole vertices)
     *
     * This method uses an O(#edges) explicit computation and should only be used for checking
     * validity
     *
     * @return number of vertices in the mesh
     */
    long count_vertices_from_orbits() const;

    /**
     * @brief Count the number of face in the mesh (including hole faces)
     *
     * This method uses an O(#edges) explicit computation and should only be used for checking
     * validity
     *
     * @return number of faces in the mesh
     */
    long count_faces_from_orbits() const;
};

} // namespace wmtk::tests
