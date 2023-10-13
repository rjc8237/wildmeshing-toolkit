#include "polygon_mesh_topology_initialization.h"
#include <Eigen/Sparse>

namespace wmtk {
namespace {

// Flatten a list of lists to a single vector
VectorXl flatten(const std::vector<std::vector<long>>& list_of_lists)
{
    // Get size of flattened array
    long arr_size = 0;
    for (long i = 0; i < list_of_lists.size(); ++i) {
        arr_size += list_of_lists[i].size();
    }

    // Build flattened array
    VectorXl arr(arr_size);
    long count = 0;
    for (long i = 0; i < list_of_lists.size(); ++i) {
        for (long j = 0; j < list_of_lists[i].size(); ++j) {
            arr[count] = list_of_lists[i][j];
            count++;
        }
    }
    return arr;
}

// Generate a range vector of integers [start, end)
std::vector<long> range(long start, long end)
{
    assert(start <= end);

    std::vector<long> arr(end - start);
    for (long i = 0; i < end - start; ++i) {
        arr[i] = start + i;
    }

    return arr;
}

// Build a right inverse g for a function f so that f(g(i)) = i
VectorXl build_right_inverse(Eigen::Ref<const VectorXl> f)
{
    long n = 0;
    for (long i = 0; i < f.size(); ++i) {
        n = std::max<long>(n, f[i] + 1);
    }
    VectorXl g(n);
    for (long i = 0; i < f.size(); ++i) {
        g[f[i]] = i;
    }

    return g;
}

// Build an inverse g for a function f
VectorXl build_inverse(Eigen::Ref<const VectorXl> f)
{
    long n = f.size();
    VectorXl g(n);
    for (long i = 0; i < n; ++i) {
        g[f[i]] = i;
    }

    return g;
}

// Permute an index map m on the left by a permutation g (i.e., i -> g[m[i]])
VectorXl reindex_map_image(Eigen::Ref<const VectorXl> map, Eigen::Ref<const VectorXl> perm)
{
    long domain_size = map.size();
    VectorXl permuted_map(domain_size);
    for (long i = 0; i < domain_size; ++i) {
        permuted_map[i] = perm[map[i]];
    }

    return permuted_map;
}

// Permute an index map m on the right by a permutation g (i.e., i -> m[g^{-1}[i]]))
VectorXl reindex_map_domain(Eigen::Ref<const VectorXl> map, Eigen::Ref<const VectorXl> perm)
{
    long domain_size = map.size();
    VectorXl permuted_map(domain_size);
    for (long i = 0; i < domain_size; ++i) {
        permuted_map[perm[i]] = i;
    }

    return permuted_map;
}

// Conjugate an index map m on the by a permutation g (i.e., i -> g[m[g^{-1}[i]]]))
VectorXl conjugate_permutation_map(Eigen::Ref<const VectorXl> map, Eigen::Ref<const VectorXl> perm)
{
    long domain_size = map.size();
    VectorXl permuted_map(domain_size);
    for (long i = 0; i < domain_size; ++i) {
        permuted_map[perm[i]] = perm[map[i]];
    }

    return permuted_map;
}

// Build the opposite halfedge map from halfedge to tail (from) and head (to) vertex maps,
// gluing together halfedges with common endpoints
VectorXl build_opp(Eigen::Ref<const VectorXl> from, Eigen::Ref<const VectorXl> to)
{
    assert(to.size() == from.size());
    long n_he = from.size();
    long n_v = 0;
    for (long i = 0; i < from.size(); ++i) {
        n_v = std::max<long>(n_v, from[i] + 1);
    }

    // Shift indices by 1 to distinguish from 0 entries in the matrix
    std::vector<long> he_index = range(1, n_he + 1);

    // Create the adjacency matrix (from, to) -> halfedge index+1
    Eigen::SparseMatrix<long> vv2he(n_v, n_v);
    typedef Eigen::Triplet<long> Trip;
    std::vector<Trip> trips(n_he);
    for (long he = 0; he < n_he; ++he) {
        trips[he] = Trip(from[he], to[he], he_index[he]);
    }
    vv2he.setFromTriplets(trips.begin(), trips.end());

    // Create opp array
    // TODO Check validity and manifold
    VectorXl opp(n_he);
    for (long he = 0; he < n_he; ++he) {
        opp[he] = vv2he.coeffRef(to[he], from[he]) - 1;
    }

    return opp;
}

// Extend next and opp to add extra halfedges along the boundaries.
std::tuple<VectorXl, VectorXl, VectorXl, VectorXl, std::vector<long>> build_boundary_loops(
    Eigen::Ref<const VectorXl> next,
    Eigen::Ref<const VectorXl> opp,
    Eigen::Ref<const VectorXl> to,
    Eigen::Ref<const VectorXl> he2f)
{
    long n_he = opp.size();

    // Get halfedges on the boundary
    std::vector<long> bnd_he;
    for (long he = 0; he < n_he; ++he) {
        if (opp[he] == -1) {
            bnd_he.push_back(he);
        }
    }
    long n_bnd_he = bnd_he.size();

    // Add boundary halfedges to the extended opp array
    VectorXl opp_ext(n_he + n_bnd_he);
    for (long i = 0; i < n_he; ++i) {
        opp_ext[i] = opp[i];
    }
    for (long i = 0; i < n_bnd_he; ++i) {
        opp_ext[bnd_he[i]] = n_he + i;
        opp_ext[n_he + i] = bnd_he[i];
    }

    // Add boundary halfedge data for the extended next and to halfedge array
    VectorXl next_ext(n_he + n_bnd_he);
    VectorXl to_ext(n_he + n_bnd_he);
    for (long i = 0; i < n_he; ++i) {
        next_ext[i] = next[i];
        to_ext[i] = to[i];
    }
    for (long i = 0; i < n_bnd_he; ++i) {
        long he = bnd_he[i];
        long he_it = next[he];
        while (opp[he_it] != -1) {
            he_it = next[opp[he_it]];
        }
        next_ext[opp_ext[he_it]] = opp_ext[he];
        to_ext[opp_ext[he_it]] = to[he];
    }

    // Build faces for extended halfedge, including new boundary faces
    std::vector<std::vector<long>> faces = build_orbits(next_ext);
    long n_f_ext = faces.size();
    long n_f = 0;
    for (long i = 0; i < he2f.size(); ++i) {
        n_f = std::max<long>(n_f, he2f[i] + 1);
    }

    // Add boundary halfedge data for the extended he2f array, recording a halfedge per new face
    std::vector<long> bnd_loops(0);
    VectorXl he2f_ext(n_he + n_bnd_he);
    for (long i = 0; i < n_he; ++i) {
        he2f_ext[i] = he2f[i];
    }
    for (long i = 0; i < n_f_ext; ++i) {
        if (faces[i][0] >= n_he) {
            bnd_loops.push_back(faces[i][0]);
            long face_size = faces[i].size();
            for (long j = 0; j < face_size; ++j) {
                he2f_ext[faces[i][j]] = n_f;
            }
            ++n_f; // Increment number of faces
        }
    }

    return std::make_tuple(next_ext, opp_ext, to_ext, he2f_ext, bnd_loops);
}

// Generate halfedge permutation to make opp implicit
VectorXl build_implicit_opp_reindex(Eigen::Ref<const VectorXl> opp)
{
    long n_he = opp.size();

    // Get edge pairs
    std::vector<std::vector<long>> edges = build_orbits(opp);
    long n_e = edges.size();
    assert(2 * n_e == n_he);

    // Build permutation mapping old halfedge indices to new indices with implicit opp
    VectorXl he_reindex(n_he);
    for (long e = 0; e < n_e; ++e) {
        assert(edges[e].size() == 2);
        he_reindex[edges[e][0]] = 2 * e;
        he_reindex[edges[e][1]] = 2 * e + 1;
    }

    return he_reindex;
}

// Implicit opposite map for halfedges paired as [e] = {2*e, 2*e + 1}
long implicit_opp(long h)
{
    return ((h % 2) == 0) ? (h + 1) : (h - 1);
}

} // namespace

std::vector<std::vector<long>> build_orbits(Eigen::Ref<const VectorXl> perm)
{
    long n_perm = perm.size();

    // Get the maximum value in perm
    long max_perm = 0;
    for (long i = 0; i < n_perm; ++i) {
        max_perm = std::max(perm[i], max_perm);
    }
    std::vector<bool> visited(max_perm + 1, false);

    // Build the cycles of the permutation
    std::vector<std::vector<long>> cycles(0);
    for (long i = 0; i < max_perm + 1; ++i) {
        if (!visited[i]) {
            cycles.push_back(std::vector<long>());
            long i_it = i;
            while (true) {
                visited[i_it] = true;
                cycles.back().push_back(i_it);
                i_it = perm[i_it];
                if (i_it == i) {
                    break;
                }
            }
        }
    }

    return cycles;
}

std::tuple<VectorXl, VectorXl, VectorXl, std::vector<long>> fv_to_nh(
    std::vector<std::vector<long>>& F)
{
    // Get the cumulative sum of the number of halfedges per face
    long n_f = F.size();
    long n_he = 0;
    std::vector<long> F_he_cumsum(n_f);
    for (long i = 0; i < n_f; ++i) {
        n_he += F[i].size();
        F_he_cumsum[i] = n_he;
    }

    // Create a list of list of indices of halfedges per face, not including boundary-loop faces.
    // Halfedges of face [v_1,...,v_n] are ordered (v_1 -> v_2),...,(v_{n-1} -> v_n),(v_n -> v_1)
    std::vector<std::vector<long>> F_he(n_f);
    F_he[0] = range(0, F_he_cumsum[0]);
    for (long i = 1; i < n_f; ++i) {
        F_he[i] = range(F_he_cumsum[i - 1], F_he_cumsum[i]);
    }

    // Create the per face next halfedge map
    std::vector<std::vector<long>> F_n(n_f);
    for (long i = 0; i < n_f; ++i) {
        F_n[i] = std::vector<long>(F_he[i].size());
        for (long j = 0; j < F_he[i].size() - 1; ++j) {
            F_n[i][j] = F_he[i][j + 1];
        }
        F_n[i][F_he[i].size() - 1] = F_he[i][0];
    }

    // Get the per face indices of head vertices of halfedges
    std::vector<std::vector<long>> F_to(n_f);
    for (long i = 0; i < n_f; ++i) {
        F_to[i] = std::vector<long>(F[i].size());
        for (long j = 0; j < F[i].size() - 1; ++j) {
            F_to[i][j] = F[i][j + 1];
        }
        F_to[i][F[i].size() - 1] = F[i][0];
    }

    // Get the per face indices of adjacent faces of halfedges
    std::vector<std::vector<long>> F_face(n_f);
    for (long i = 0; i < n_f; ++i) {
        F_face[i] = std::vector<long>(F[i].size(), i);
    }

    // Flatten the per face maps to per halfedge maps
    // WARNING: Assumes the flatten operation transforms F_he to the identity halfedge map
    VectorXl next = flatten(F_n);
    VectorXl from = flatten(F); // the tail of the halfedge (f, hi) is just vi
    VectorXl to = flatten(F_to);
    VectorXl he2f = flatten(F_face);

    // Generate opposite halfedge map
    VectorXl opp = build_opp(from, to);

    // Create boundary loops
    std::vector<long> bnd_loops;
    std::tie(next, opp, to, he2f, bnd_loops) = build_boundary_loops(next, opp, to, he2f);
    long n_bd_f = bnd_loops.size();

    // Reindex halfedges to make opp implicit, i.e., make adjacent halfedges opposite
    VectorXl he_reindex = build_implicit_opp_reindex(opp);
    conjugate_permutation_map(next, he_reindex);
    reindex_map_domain(to, he_reindex);
    reindex_map_domain(he2f, he_reindex);
    for (long i = 0; i < n_bd_f; ++i) {
        bnd_loops[i] = he_reindex[bnd_loops[i]];
    }

    return std::make_tuple(next, to, he2f, bnd_loops);
}

std::tuple<VectorXl, VectorXl, VectorXl, std::vector<long>> fv_to_nh(
    Eigen::Ref<const RowVectors3l> F)
{
    // Convert eigen matrix to vector of vectors
    std::vector<std::vector<long>> F_vec(F.rows(), std::vector<long>(F.cols()));
    for (long i = 0; i < F.rows(); ++i) {
        for (long j = 0; j < F.cols(); ++j) {
            F_vec[i][j] = F(i, j);
        }
    }

    // Run fh_to_nh method
    return fv_to_nh(F_vec);
}

std::tuple<VectorXl, VectorXl, VectorXl, VectorXl, VectorXl> polygon_mesh_topology_initialization(
    Eigen::Ref<const VectorXl> next)
{
    // Build previous halfedge array
    long n_he = next.size();
    VectorXl prev = build_inverse(next);

    // Construct face loops
    std::vector<std::vector<long>> faces = build_orbits(next);

    // Map faces to attached halfedges
    long n_f = faces.size();
    VectorXl f2he(n_f);
    for (long i = 0; i < n_f; ++i) {
        f2he[i] = faces[i][0];
    }

    // Map halfedges to faces
    std::vector<long> he2f_perm;
    he2f_perm.reserve(n_he);
    for (long i = 0; i < n_f; ++i) {
        for (long j = 0; j < faces[i].size(); ++j) {
            he2f_perm.push_back(i);
        }
    }
    VectorXl Fhe = flatten(faces);
    VectorXl he2f(n_he);
    for (long i = 0; i < n_he; ++i) {
        he2f[Fhe[i]] = he2f_perm[i];
    }

    // Create vertex list
    VectorXl circ(next.size());
    for (long he = 0; he < next.size(); ++he) {
        circ[he] = prev[implicit_opp(he)];
    }
    std::vector<std::vector<long>> vert = build_orbits(circ);

    // Create out array
    long n_v = vert.size();
    VectorXl out(n_v);
    for (long i = 0; i < n_v; ++i) {
        out[i] = next[vert[i][0]];
    }

    // Create to array
    std::vector<long> vind;
    vind.reserve(n_he);
    for (long i = 0; i < n_v; ++i) {
        for (long j = 0; j < vert[i].size(); ++j) {
            vind.push_back(i);
        }
    }
    VectorXl vhe = flatten(vert);
    VectorXl to(n_he);
    for (long i = 0; i < n_he; ++i) {
        to[vhe[i]] = vind[i];
    }

    return std::make_tuple(prev, to, he2f, f2he, out);
}

std::tuple<VectorXl, VectorXl, VectorXl, VectorXl, VectorXl, VectorXl, std::vector<long>>
polygon_mesh_topology_initialization(Eigen::Ref<const RowVectors3l> F)
{
    auto [next, to, he2f, hole_faces] = fv_to_nh(F);
    VectorXl prev = build_inverse(next);

    // Build the vertex and face to halfedge array
    VectorXl f2he = build_right_inverse(he2f);
    VectorXl in = build_right_inverse(to);
    VectorXl out = reindex_map_image(in, next); // next[in[v]] emanates from v

    return std::make_tuple(next, prev, to, he2f, f2he, out, hole_faces);
}

// Deprecated: Unecessary with explicit to and he2f methods
std::tuple<std::vector<long>, std::vector<long>> build_halfedge_reindex_maps(
    Eigen::Ref<const VectorXl> next,
    Eigen::Ref<const VectorXl> to,
    Eigen::Ref<const VectorXl> he2f)
{
    // Get the number of faces and vertices
    long n_v = 0;
    for (long i = 0; i < to.size(); ++i) {
        n_v = std::max<long>(n_v, to[i] + 1);
    }
    long n_f = 0;
    for (long i = 0; i < he2f.size(); ++i) {
        n_f = std::max<long>(n_f, he2f[i] + 1);
    }

    // Create previous and vertex circulator halfedge arrays
    VectorXl prev = build_inverse(next);
    VectorXl circ(next.size());
    for (long he = 0; he < next.size(); ++he) {
        circ[he] = prev[implicit_opp(he)];
    }

    // Create vertex and face arrays from the circulator array
    std::vector<std::vector<long>> vert = build_orbits(circ);
    std::vector<std::vector<long>> face = build_orbits(next);

    // Map circ orbit indices to indices of input vertices
    std::vector<long> vtx_reindex(n_v);
    for (long i = 0; i < n_v; ++i) {
        vtx_reindex[i] = to[vert[i][0]];
    }

    // Map next orbit indices to indices of input faces
    std::vector<long> face_reindex(n_f);
    for (long i = 0; i < n_f; ++i) {
        face_reindex[i] = he2f[face[i][0]];
    }

    return std::make_tuple(vtx_reindex, face_reindex);
}


} // namespace wmtk