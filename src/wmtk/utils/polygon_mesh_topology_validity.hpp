#pragma once

#include <wmtk/Types.hpp>

namespace wmtk {

/**
 * @brief Check if a map is invariant under some permutation and thus descends to a well-defined
 * function on the orbits of the permutation
 *
 * @param map: map from {0,...,n-1} to {0,...,m-1}
 * @param perm: permutation of n elements
 * @return true iff the map is invariant under perm
 */
bool is_invariant_under_permutation(
    Eigen::Ref<const VectorXl> map,
    Eigen::Ref<const VectorXl> perm);

/**
 * @brief Check if two maps are one-sided inverses of each other
 *
 * @param left_inverse: map from {0,...,m-1} to {0,...,n-1}
 * @param right_inverse: map from {0,...,n-1} to {0,...,m-1}
 * @return true iff left_inverse composed with right_inverse is the identity
 */
bool is_one_sided_inverse(
    Eigen::Ref<const VectorXl> left_inverse,
    Eigen::Ref<const VectorXl> right_inverse);


/**
 * @brief Check if the maps defining the edge connectivity of a polygonal mesh (next and prev)
 * are valid
 *
 * @param next: size #he vector, next halfedge id
 * @param prev: size #he vector, prev halfedge id
 * @return true iff the edge maps are valid
 */
bool are_polygon_mesh_edges_valid(Eigen::Ref<const VectorXl> next, Eigen::Ref<const VectorXl> prev);


/**
 * @brief Check if the maps defining the vertex connectivity of a polygonal mesh (to and out)
 * are valid
 *
 * @param prev: size #he vector, prev halfedge id
 * @param to: size #he vector, halfedge vertex tip id
 * @param out: size #v vector, arbitrary halfedge id outgoing from vertex
 * @return true iff the vertex maps are valid
 */
bool are_polygon_mesh_vertices_valid(
    Eigen::Ref<const VectorXl> prev,
    Eigen::Ref<const VectorXl> to,
    Eigen::Ref<const VectorXl> out);

/**
 * @brief Check if the maps defining the face connectivity of a polygonal mesh (he2f and f2he)
 * are valid
 *
 * @param next: size #he vector, next halfedge id
 * @param he2f: size #he vector, face id adjacent to halfedge
 * @param f2he: size #f vector, arbitrary halfedge id adjacent to face
 * @return true iff the face maps are valid
 */
bool are_polygon_mesh_faces_valid(
    Eigen::Ref<const VectorXl> next,
    Eigen::Ref<const VectorXl> he2f,
    Eigen::Ref<const VectorXl> f2he);

} // namespace wmtk