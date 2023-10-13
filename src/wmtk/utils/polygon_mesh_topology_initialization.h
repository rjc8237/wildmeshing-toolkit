
#pragma once

#include <wmtk/Types.hpp>


namespace wmtk {

/**
 * @brief Build orbits following next id recorded in perm.
 *
 * @param perm: a permutation given by a list of non-repeating integers in the range 0..len(perm)
 * @return cycles: a list of lists, each list represents a cycle of perm
 */
std::vector<std::vector<long>> build_orbits(Eigen::Ref<const VectorXl> perm);

/**
 * @brief Convert from list of faces (FV) mesh representation to NH data structure
 *
 * Opp is implicit with halfedges for edge e given by 2 * e and 2 * e + 1
 *
 * @param F: each row represents the vertex id (ccw) of current face
 * @return next: (N) size #h vector, next halfedge id
 * @return to: size #h vector, tip vertex of halfedge id
 * @return he2f: size #h vector, adjacent face of halfedge id
 * @return hole_faces: (H) collection of boundary face halfedge ids.
 */
std::tuple<VectorXl, VectorXl, VectorXl, std::vector<long>> fv_to_nh(
    std::vector<std::vector<long>>& F);

/**
 * @brief Convert from a FV triangle mesh to NH data structure
 */
std::tuple<VectorXl, VectorXl, VectorXl, std::vector<long>> fv_to_nh(
    Eigen::Ref<const RowVectors3l> F);

/**
 * @brief Convert from minimal next (and implicit opp) halfedge representation to full halfedge
 * surface connectivity
 *
 * @param next: size #he vector, next halfedge id
 * @return prev: size #he vector, prev halfedge id
 * @return to: size #he vector, halfedge vertex tip id
 * @return out: size #v vector, arbitrary halfedge id outgoing from vertex
 * @return he2f: size #he vector, face id adjacent to halfedge
 * @return f2he: size #f vector, arbitrary halfedge id adjacent to face
 */
std::tuple<VectorXl, VectorXl, VectorXl, VectorXl, VectorXl> polygon_mesh_topology_initialization(
    Eigen::Ref<const VectorXl> next);

/**
 * @brief Convert from matrix (F) mesh representation to full halfedge connectivity
 *
 * @param F: #f*3, each row represents the vertex id (ccw) of current face
 * @return next: size #he vector, next halfedge id
 * @return prev: size #he vector, prev halfedge id
 * @return to: size #he vector, halfedge vertex tip id
 * @return out: size #v vector, arbitrary halfedge id outgoing from vertex
 * @return he2f: size #he vector, face id adjacent to halfedge
 * @return f2he: size #f vector, arbitrary halfedge id adjacent to face
 * @return hole_faces: collection of boundary face halfedge ids.
 */
std::tuple<VectorXl, VectorXl, VectorXl, VectorXl, VectorXl, VectorXl, std::vector<long>>
polygon_mesh_topology_initialization(Eigen::Ref<const RowVectors3l> F);


} // namespace wmtk