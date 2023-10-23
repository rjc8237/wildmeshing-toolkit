#include "PolygonMesh_examples.hpp"

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
    F[0] = std::vector<long>({0, 1, 2, 3});
    F[1] = std::vector<long>({3, 2, 4});
    m.initialize_fv(F);
    return m;
}

PolygonMesh grid()
{
    PolygonMesh m;
    assert(false);
    // TODO Implement with FE representation
    return m;
}


} // namespace wmtk::tests
