#include "polygon_mesh_topology_validity.hpp"

#include "polygon_mesh_topology_initialization.h"

namespace wmtk {

namespace {
// Implicit opposite map for halfedges paired as [e] = {2*e, 2*e + 1}
long implicit_opp(long h)
{
    return ((h % 2) == 0) ? (h + 1) : (h - 1);
}
} // namespace


bool is_invariant_under_permutation(Eigen::Ref<const VectorXl> map, Eigen::Ref<const VectorXl> perm)
{
    if (map.size() != perm.size()) {
        return false;
    }
    long n = map.size();
    for (long i = 0; i < n; ++i) {
        if ((perm(i) < 0) || (perm(i) >= n) || map(perm(i)) != map(i)) {
            return false;
        }
    }
    return true;
}


bool is_one_sided_inverse(
    Eigen::Ref<const VectorXl> left_inverse,
    Eigen::Ref<const VectorXl> right_inverse)
{
    long n = right_inverse.size();
    long m = left_inverse.size();
    for (long i = 0; i < n; ++i) {
        if ((right_inverse(i) < 0) || (right_inverse(i) >= m) ||
            (left_inverse(right_inverse(i)) != i)) {
            return false;
        }
    }
    return true;
}


bool are_polygon_mesh_edges_valid(Eigen::Ref<const VectorXl> next, Eigen::Ref<const VectorXl> prev)
{
    if (next.size() != prev.size()) {
        return false;
    }

    // prev is a right and left inverse for next
    if ((!is_one_sided_inverse(next, prev)) || (!is_one_sided_inverse(prev, next))) {
        return false;
    }

    return true;
}


bool are_polygon_mesh_vertices_valid(
    Eigen::Ref<const VectorXl> prev,
    Eigen::Ref<const VectorXl> to,
    Eigen::Ref<const VectorXl> out)
{
    if (prev.size() != to.size()) {
        return false;
    }
    long n_halfedges = to.size();

    // Generate per halfedge vertex circulation and from maps
    VectorXl circ(n_halfedges);
    VectorXl from(n_halfedges);
    for (long hi = 0; hi < n_halfedges; ++hi) {
        circ(hi) = prev(implicit_opp(hi));
        from(hi) = to(implicit_opp(hi));
    }

    // Build vertices from vertex circulation
    std::vector<std::vector<long>> vert = build_orbits(circ);
    long n_vertices = vert.size();

    // Number of vertices in out match the number of orbits
    if (out.size() != n_vertices) {
        return false;
    }

    // to is invariant under circulation
    if (!is_invariant_under_permutation(to, circ)) {
        return false;
    }

    // out is a right inverse for from
    if (!is_one_sided_inverse(from, out)) {
        return false;
    }

    return true;
}

bool are_polygon_mesh_faces_valid(
    Eigen::Ref<const VectorXl> next,
    Eigen::Ref<const VectorXl> he2f,
    Eigen::Ref<const VectorXl> f2he)
{
    if (next.size() != he2f.size()) {
        return false;
    }
    long n_halfedges = he2f.size();

    // Build faces from next map
    std::vector<std::vector<long>> faces = build_orbits(next);
    long n_faces = faces.size();

    // Get number of boundary faces (i.e., number of negative indices)
    long n_bd_faces = std::abs(he2f.minCoeff());

    // Number of faces in f2he and boundary faces match the number of orbits
    if ((f2he.size() + n_bd_faces) != n_faces) {
        return false;
    }

    // he2f is invariant under next
    if (!is_invariant_under_permutation(he2f, next)) {
        return false;
    }

    // f2he is a right inverse for he2f
    if (!is_one_sided_inverse(he2f, f2he)) {
        return false;
    }

    return true;
}

} // namespace wmtk