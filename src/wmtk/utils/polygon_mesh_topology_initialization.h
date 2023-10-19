
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
 * @brief Convert from face vertex list (F) mesh representation to full halfedge connectivity
 *
 * Opp is implicit with halfedges for edge e given by 2 * e and 2 * e + 1
 *
 * @param F: each row represents the vertex id (ccw) of current face
 * @return next: size #he vector, next halfedge id
 * @return prev: size #he vector, prev halfedge id
 * @return to: size #he vector, halfedge vertex tip id
 * @return out: size #v vector, arbitrary halfedge id outgoing from vertex
 * @return he2f: size #he vector, face id adjacent to halfedge
 * @return f2he: size #f vector, arbitrary halfedge id adjacent to face
 * @return hole_faces: collection of boundary face halfedge ids.
 */
std::tuple<VectorXl, VectorXl, VectorXl, VectorXl, VectorXl, VectorXl, std::vector<long>>
polygon_mesh_topology_initialization(std::vector<std::vector<long>>& F);

std::tuple<VectorXl, VectorXl, VectorXl, VectorXl, VectorXl, VectorXl, std::vector<long>>
polygon_mesh_topology_initialization(Eigen::Ref<const RowVectors3l> F);


} // namespace wmtk