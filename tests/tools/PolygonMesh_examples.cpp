#include "PolygonMesh_examples.hpp"

namespace wmtk::tests {


//        * \ 
//       /   e1
//      e2    \ 
//     /       *
//     *---e0--->
//
PolygonMesh triangle()
{
    PolygonMesh m;

    VectorXl next(6);
    next << 2, 5, 4, 1, 0, 3;
    std::vector<long> bnd_loops = {1};
    m.initialize(next, bnd_loops);
    return m;
}

//     <---e2---*
//    *          |
//    |          e1
//    e3         |
//    |          *
//     *---e0--->
//
PolygonMesh square()
{
    PolygonMesh m;

    VectorXl next(8);
    next << 2, 7, 4, 1, 6, 3, 0, 5;
    std::vector<long> bnd_loops = {1};
    m.initialize(next, bnd_loops);
    return m;
}

//          * \ 
//         /    e2
//      e3         \ 
//     /            *
//   *               /
//    \             e1
//     e4          /
//       \        *
//        *--e0-->
//
PolygonMesh pentagon()
{
    PolygonMesh m;

    VectorXl next(10);
    next << 2, 9, 4, 1, 6, 3, 8, 5, 0, 7;
    std::vector<long> bnd_loops = {1};
    m.initialize(next, bnd_loops);
    return m;
}

//       -*<-e0-
//     /         \ 
//    |           |
//    |           |
//     \         /
//       -e1->*-
//
PolygonMesh digon()
{
    PolygonMesh m;

    VectorXl next(4);
    next << 2, 3, 0, 1;
    std::vector<long> bnd_loops = {1};
    m.initialize(next, bnd_loops);
    return m;
}

//      --e0--
//    /        \ 
//   |          |
//   |          |
//    \        /
//      -->*--
//
PolygonMesh monogon()
{
    PolygonMesh m;

    VectorXl next(2);
    next << 0, 1;
    std::vector<long> bnd_loops = {1};
    m.initialize(next, bnd_loops);
    return m;
}

//      -------
//    /         \ 
//   |  *-h0a->  |
//   |  <-h0b-*  |
//    \         /
//      -------
PolygonMesh bubble()
{
    PolygonMesh m;

    VectorXl next(2);
    next << 1, 0;
    std::vector<long> bnd_loops = {};
    m.initialize(next, bnd_loops);
    return m;
}

//      -------
//    /         \ 
//   |-->*--h0a--|
//   |--*<--h0b--|
//    \         /
//      -------
PolygonMesh dual_bubble()
{
    PolygonMesh m;

    VectorXl next(2);
    next << 0, 1;
    std::vector<long> bnd_loops = {};
    m.initialize(next, bnd_loops);
    return m;
}

//     <---e8---*    <---e11--*
//    *          |  *          |
//    |         h7a |          e10
//    e9         | h7b         |
//    |          *  |          *
//     *---h2b-->    *---e6b-->
//     <---h2a--*    <---h6a--*
//    *          |  *          |
//    |         h1a |          e5
//    e3         | h1b         |
//    |          *  |          *
//     *---e0--->    *---e4--->
//
PolygonMesh grid()
{
    PolygonMesh m;
    assert(false);
    // TODO Implement with FE representation
    return m;
}


} // namespace wmtk::tests
