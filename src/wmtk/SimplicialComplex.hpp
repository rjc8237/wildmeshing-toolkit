#pragma once

#include <cassert>
#include <queue>
#include <set>
#include <vector>
#include "Simplex.hpp"
#include "TetMesh.hpp"
#include "TriMesh.hpp"

using namespace wmtk;

class SimplicialComplex
{
private:
    std::set<Simplex> simplexes;

public:
    const std::set<Simplex>& get_simplices() const { return simplexes; }

    std::set<Simplex> get_simplices(const PrimitiveType& ptype) const
    {
        std::set<Simplex> ret;
        for (const Simplex& s : simplexes) {
            if (s.primitive_type() == ptype) {
                ret.insert(s);
            }
        }

        return ret;
    }

    /**
     * @brief Add simplex to the complex if it is not already in it.
     *
     * @returns false if simplex is already in the complex
     */
    bool add_simplex(const Simplex& s)
    {
        assert(s.primitive_type() != PrimitiveType::Invalid);
        const auto [it, was_successful] = simplexes.insert(s);
        return was_successful;
    }

    void unify_with_complex(const SimplicialComplex& other)
    {
        // this is N log(N) complexity
        for (const Simplex& s : other.get_simplices()) {
            add_simplex(s);
        }
    }

    bool operator==(const SimplicialComplex& other) const
    {
        if (simplexes.size() != other.simplexes.size()) {
            return false;
        }
        // this is N log(N) complexity
        for (const auto& t1 : simplexes) {
            const auto it = other.simplexes.find(t1);
            if (it == other.simplexes.end()) {
                return false;
            }
        }

        return true;
    }

    SimplicialComplex& operator=(const SimplicialComplex&) = default;

    SimplicialComplex() = default;

    SimplicialComplex(const std::vector<Simplex>& ss)
    {
        for (const Simplex& s : ss) {
            add_simplex(s);
        }
    }
};

inline SimplicialComplex get_union(const SimplicialComplex& sc1, const SimplicialComplex& sc2)
{
    SimplicialComplex u = sc1;
    u.unify_with_complex(sc2);
    return u;
}

inline SimplicialComplex get_intersection(const SimplicialComplex& A, const SimplicialComplex& B)
{
    SimplicialComplex sc_union = A;
    SimplicialComplex sc_intersection;

    for (const auto& s : B.get_simplices()) {
        if (!sc_union.add_simplex(s)) {
            // s is already in A --> s is in the intersection of A and B
            sc_intersection.add_simplex(s);
        }
    }

    return sc_intersection;
}

//////////////////////////////////
// List of Operators
// bd: boundary
// clbd: closed boudnary
// st: start
// clst: closed star
// lnk: link
//////////////////////////////////

// ∂s
/**
 * @brief get the boundary of a simplex
 */
SimplicialComplex boundary(const Simplex& s, const Mesh& m)
{
    SimplicialComplex sc;

    constexpr PrimitiveType PV = PrimitiveType::Vertex;
    constexpr PrimitiveType PE = PrimitiveType::Edge;
    constexpr PrimitiveType PF = PrimitiveType::Face;

    const Tuple t = s.tuple();

    auto sw = [&m](const Tuple& t, const PrimitiveType& ptype) { return m.switch_tuple(t, ptype); };


    // exhaustive implementation
    switch (s.primitive_type()) {
    case PrimitiveType::Tetrahedron:
        // bd(tet) = 4triangles + 6 edges + 4vertices
        sc.add_simplex(Simplex(PV, t)); // A
        sc.add_simplex(Simplex(PV, sw(t, PV))); // B
        sc.add_simplex(Simplex(PV, sw(sw(t, PE),
                                      PV))); // C
        sc.add_simplex(Simplex(PV, sw(sw(t, PF),
                                      PV))); // D
        sc.add_simplex(Simplex(PE, t)); // AB
        sc.add_simplex(Simplex(PE, sw(t, PE))); // AC
        sc.add_simplex(Simplex(PE, sw(sw(t, PV),
                                      PE))); // BC
        sc.add_simplex(Simplex(PE, sw(sw(t, PF),
                                      PE))); // AD
        sc.add_simplex(Simplex(PE, sw(sw(sw(t, PV), PF),
                                      PE))); // BD
        sc.add_simplex(Simplex(PE, sw(sw(sw(sw(t, PE), PV), PF),
                                      PE))); // CD
        sc.add_simplex(Simplex(PF, t)); // ABC
        sc.add_simplex(Simplex(PF, sw(t, PF))); // ABD
        sc.add_simplex(Simplex(PF, sw(sw(t, PE),
                                      PF))); // ACD
        sc.add_simplex(Simplex(PF, sw(sw(sw(t, PV), PE),
                                      PF))); // BCD
        break;
    case PF: // bd(triangle) = 3edges + 3vertices
        sc.add_simplex(Simplex(PV, t));
        sc.add_simplex(Simplex(PV, sw(t, PV)));
        sc.add_simplex(Simplex(PV, sw(sw(t, PE), PV)));
        sc.add_simplex(Simplex(PE, t));
        sc.add_simplex(Simplex(PE, sw(t, PE)));
        sc.add_simplex(Simplex(PE, sw(sw(t, PV), PE)));
        /* code */
        break;
    case PE:
        // bd(edge) = 2 vertices
        sc.add_simplex(Simplex(PV, t));
        sc.add_simplex(Simplex(PV, sw(t, PV)));
        /* code */
        break;
    case PV: break;
    default: assert(false); break;
    }

    return sc;
}

// ∂s∪{s}
/**
 * @brief get complex of a simplex and its boundary
 */
SimplicialComplex simplex_with_boundary(const Simplex& s, const Mesh& m)
{
    SimplicialComplex sc = boundary(s, m);
    sc.add_simplex(s);
    return sc;
}

// Simplex s1,s2, check if A∩B!=∅
// check is intersect(∂s1, ∂s2) has intersections
/**
 * @brief check if simplices with their boundary intersect
 */
inline bool simplices_w_boundary_intersect(const Simplex& s1, const Simplex& s2, const Mesh& m)
{
    SimplicialComplex s1_bd = simplex_with_boundary(s1, m);
    SimplicialComplex s2_bd = simplex_with_boundary(s2, m);
    SimplicialComplex s1_s2_int = get_intersection(s1_bd, s2_bd);
    return (s1_s2_int.get_simplices().size() != 0);
}

SimplicialComplex closed_star(const Simplex& s, const Mesh& m)
{
    SimplicialComplex sc;

    constexpr PrimitiveType PV = PrimitiveType::Vertex;
    constexpr PrimitiveType PE = PrimitiveType::Edge;
    constexpr PrimitiveType PF = PrimitiveType::Face;
    constexpr PrimitiveType PT = PrimitiveType::Tetrahedron;

    const Tuple t = s.tuple();

    auto sw = [&m](const Tuple& t, const PrimitiveType& ptype) { return m.switch_tuple(t, ptype); };

    // const int &cell_dim = m->cell_dimension(); // TODO: 2 for trimesh, 3 for tetmesh need it in Mesh class
    const int cell_dim = dynamic_cast<const TriMesh*>(&m) ? 2 : 3;
    if (cell_dim == 2) {
        switch (s.primitive_type()) {
        case PV: {
            std::queue<Tuple> q;
            q.push(s.tuple());
            while (!q.empty()) {
                const Tuple t = q.front();
                q.pop();
                if (sc.add_simplex(Simplex(PF, t))) {
                    if (!m.is_boundary(t)) {
                        q.push(sw(t, PF));
                    }
                    if (!m.is_boundary(sw(t, PE))) {
                        q.push(sw(sw(t, PE), PF));
                    }
                }
            }
            break;
        }
        case PE:
            sc.add_simplex(Simplex(PF, t));
            if (!m.is_boundary(t)) {
                sc.add_simplex(Simplex(PF, sw(t, PF)));
            }
            break;
        case PF: sc.add_simplex(s); break;
        default: assert(false); break;
        }
    } else if (cell_dim == 3) {
        switch (s.primitive_type()) {
        case PV: {
            std::queue<Tuple> q;
            q.push(t);
            while (!q.empty()) {
                Tuple t = q.front();
                q.pop();
                if (sc.add_simplex(Simplex(PF, t))) {
                    const Tuple t1 = t;
                    const Tuple t2 = sw(t, PF);
                    const Tuple t3 = sw(sw(t, PE), PF);
                    if (!m.is_boundary(t1)) {
                        q.push(sw(t1, PT));
                    }
                    if (!m.is_boundary(t2)) {
                        q.push(sw(t2, PT));
                    }
                    if (!m.is_boundary(t3)) {
                        q.push(sw(t3, PT));
                    }
                }
            }
            break;
        }
        case PE: {
            std::queue<Tuple> q;
            q.push(t);
            while (!q.empty()) {
                Tuple t = q.front();
                q.pop();
                if (sc.add_simplex(Simplex(PT, t))) {
                    if (!m.is_boundary(t)) {
                        q.push(sw(t, PT));
                    }
                    if (!m.is_boundary(sw(t, PF))) {
                        q.push(sw(sw(t, PF), PT));
                    }
                }
            }
            break;
        }
        case PF: {
            sc.add_simplex(Simplex(PT, t));
            if (!m.is_boundary(t)) {
                sc.add_simplex(Simplex(PT, sw(t, PT)));
            }
            break;
        }
        case PT: {
            sc.add_simplex(s);
            break;
        }
        default: {
            assert(false);
            break;
        }
        }
    }

    const auto top_simplices = sc.get_simplices();
    for (const Simplex& ts : top_simplices) {
        sc.unify_with_complex(boundary(ts, m));
    }
    return sc;
}

SimplicialComplex link(const Simplex& s, const Mesh& m)
{
    SimplicialComplex sc_clst = closed_star(s, m);
    SimplicialComplex sc;
    for (const Simplex& ss : sc_clst.get_simplices()) {
        if (!simplices_w_boundary_intersect(s, ss, m)) {
            sc.add_simplex(ss);
        }
    }

    return sc;
}

SimplicialComplex open_star(const Simplex& s, const Mesh& m)
{
    SimplicialComplex sc_clst = closed_star(s, m);
    SimplicialComplex sc;
    sc.add_simplex(s);
    for (const Simplex& ss : sc_clst.get_simplices()) {
        if (ss.primitive_type() <= s.primitive_type()) {
            continue;
        }
        if (simplices_w_boundary_intersect(s, ss, m)) {
            sc.add_simplex(ss);
        }
    }

    return sc;
}

//////////////////////////////////
// check link condition
// input Tuple t --> edge (a,b)
// check if lnk(a) ∩ lnk(b) == lnk(ab)
//////////////////////////////////
bool link_cond(Tuple t, const Mesh& m)
{
    SimplicialComplex lhs = link(Simplex(PrimitiveType::Vertex, t), m); // lnk(a)
    lhs.unify_with_complex(link(
        Simplex(PrimitiveType::Vertex, m.switch_tuple(t, PrimitiveType::Vertex)),
        m)); // Union lnk(b)

    SimplicialComplex rhs = link(Simplex(PrimitiveType::Edge, t), m); // lnk(ab)
    return (lhs == rhs);
}

//////////////////////////////////
// k-ring
//////////////////////////////////
/**
 * @brief get one ring neighbors of vertex in _t_
 */
std::vector<Simplex> vertex_one_ring(Tuple t, const Mesh& m)
{
    Simplex s(PrimitiveType::Vertex, t);
    SimplicialComplex sc_link = link(s, m);
    std::set<Simplex> one_ring_simplices = sc_link.get_simplices(PrimitiveType::Vertex);
    return std::vector<Simplex>(one_ring_simplices.begin(), one_ring_simplices.end());
}

std::vector<Simplex> k_ring(Tuple t, const Mesh& m, int k)
{
    if (k < 1) return {};

    SimplicialComplex sc(vertex_one_ring(t, m));
    for (int i = 2; i <= k; ++i) {
        const auto simplices = sc.get_simplices();
        for (const Simplex& s : simplices) {
            SimplicialComplex sc_or(vertex_one_ring(s.tuple(), m));
            sc.unify_with_complex(sc_or);
        }
    }

    std::set<Simplex> k_ring_simplices = sc.get_simplices();
    return std::vector<Simplex>(k_ring_simplices.begin(), k_ring_simplices.end());
}
