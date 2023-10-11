
#pragma once

#include <wmtk/Types.hpp>


namespace wmtk {

/**
 * @brief Convert from matrix (F) mesh representation to NH data structure
 *
 * Opp is implicit with halfedges for edge e given by 2 * e and 2 * e + 1
 *
 * @param F: #f*3, each row represents the vertex id (ccw) of current face
 * @return next: (N) size #h vector, next halfedge id
 * @return hole_faces: (H) collection of boundary face ids.
 */
std::tuple<VectorXl, VectorXl> fv_to_nh(Eigen::Ref<const RowVectors3l> F);

/**
 * @brief Extend next and opp to add extra halfedges along the boundaries.
 *
 * @param next: next-halfedge map same length as opp
 * @param opp: halfedge map, for boundary halfedges -1
 * @return next_he_ext: next_halfedge map, same length as opp_ext; newly added halfedges are linked
 *                      into boundary loops
 * @return opp_ext: opp with an extra coupled halfedge added at the end for each boundary halfedge
 *                  all halfedges have a pair
 */
std::tuple<VectorXl, VectorXl> build_boundary_loops(
    Eigen::Ref<const VectorXl> next,
    Eigen::Ref<const VectorXl> opp);

/**
 * @brief Build orbits following next id recorded in perm.
 *
 * @param perm: a permutation given by a list of non-repeating integers in the range 0..len(perm)
 * @return cycles: a list of lists, each list represents a cycle of perm
 */
std::vector<std::vector<int>> build_orbits(const std::vector<int>& perm);

/**
 * @brief Convert from minimal next (and implicit opp) halfedge representation with marked boundary
 * faces to full halfedge surface connectivity with marked boundary faces
 */
std::tuple<RowVectors2l, RowVectors2l, RowVectors2l, RowVectors2l, VectorXl, VectorXl, VectorXl>
polygon_mesh_topology_initialization(
    Eigen::Ref<const VectorXl> next,
    Eigen::Ref<const VectorXl> hole_faces);

/**
 * @brief Convert from matrix (F) mesh representation to full halfedge connectivity
 *
 * @param F: #f*3, each row represents the vertex id (ccw) of current face
 * @return next: size #e*2 matrix, next halfedge id
 * @return prev: size #e*2 matrix, prev halfedge id
 * @return EV: size #e*2 matrix, halfedge vertex tip id
 * @return EF: size #e*2 matrix, face id adjacent to halfedge
 * @return FH: size #f vector, arbitrary halfedge id adjacent to face
 * @return VH: size #v vector, arbitrary halfedge id outgoing from vertex
 * @return hole_faces: collection of boundary face ids.
 */
std::tuple<RowVectors2l, RowVectors2l, RowVectors2l, RowVectors2l, VectorXl, VectorXl, VectorXl>
polygon_mesh_topology_initialization(Eigen::Ref<const RowVectors3l> F);


} // namespace wmtk