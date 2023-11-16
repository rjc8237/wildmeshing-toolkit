#include "PolygonMesh_examples.hpp"
#include <algorithm>
#include <random>

namespace wmtk::tests {


PolygonMesh triangle()
{
    PolygonMesh m;

    VectorXl next(6);
    next << 2, 5, 4, 1, 0, 3;
    m.initialize(next);
    return m;
}

PolygonMesh square()
{
    PolygonMesh m;

    VectorXl next(8);
    next << 2, 7, 4, 1, 6, 3, 0, 5;
    m.initialize(next);
    return m;
}

PolygonMesh pentagon()
{
    PolygonMesh m;

    VectorXl next(10);
    next << 2, 9, 4, 1, 6, 3, 8, 5, 0, 7;
    m.initialize(next);
    return m;
}

PolygonMesh digon()
{
    PolygonMesh m;

    VectorXl next(4);
    next << 2, 3, 0, 1;
    m.initialize(next);
    return m;
}

PolygonMesh monogon()
{
    PolygonMesh m;

    VectorXl next(2);
    next << 0, 1;
    m.initialize(next);
    return m;
}

PolygonMesh bubble()
{
    PolygonMesh m;

    VectorXl next(2);
    next << 1, 0;
    m.initialize(next);
    return m;
}

PolygonMesh dual_bubble()
{
    PolygonMesh m;

    VectorXl next(2);
    next << 0, 1;
    m.initialize(next);
    return m;
}

PolygonMesh glued_polygons()
{
    PolygonMesh m;

    std::vector<std::vector<long>> F(2);
    F[0] = std::vector<long>({0, 1, 4, 3, 2});
    F[1] = std::vector<long>({2, 3, 4, 5});
    m.initialize_fv(F);
    return m;
}

PolygonMesh grid()
{
    PolygonMesh m;

    std::vector<std::vector<long>> F(4);
    F[0] = std::vector<long>({0, 1, 4, 3});
    F[1] = std::vector<long>({1, 2, 5, 4});
    F[2] = std::vector<long>({3, 4, 7, 6});
    F[3] = std::vector<long>({4, 5, 8, 7});
    m.initialize_fv(F);
    return m;
}

PolygonMesh torus()
{
    PolygonMesh m;

    VectorXl next(4);
    next << 2, 3, 1, 0;
    m.initialize(next);
    return m;
}

PolygonMesh monogon_with_leaf()
{
    PolygonMesh m;

    VectorXl next(4);
    next << 2, 1, 3, 0;
    m.initialize(next);
    return m;
}

PolygonMesh annulus()
{
    PolygonMesh m;

    std::vector<std::vector<long>> F(1);
    F[0] = std::vector<long>({0, 1, 3, 2, 0, 4, 6, 7, 5, 4});
    m.initialize_fv(F);
    return m;
}

PolygonMesh two_squares()
{
    PolygonMesh m;

    std::vector<std::vector<long>> F(2);
    F[0] = std::vector<long>({0, 1, 3, 2});
    F[1] = std::vector<long>({4, 5, 7, 6});
    m.initialize_fv(F);
    return m;
}

PolygonMesh random_polygon_mesh(long num_edges, long rng_seed)
{
    PolygonMesh m;

    VectorXl next(2 * num_edges);
    std::iota(next.begin(), next.end(), 0);
    std::mt19937 g(rng_seed);
    std::shuffle(next.begin(), next.end(), g);

    m.initialize(next);
    return m;
}

} // namespace wmtk::tests
