#pragma once

#include <wmtk/PolygonMesh.hpp>

namespace wmtk::tests {


//        * \ .
//       /   e1
//      e2    \ .
//     /       *
//     *---e0--->
//
PolygonMesh triangle();

//     <---e2---*
//    *          |
//    |          e1
//    e3         |
//    |          *
//     *---e0--->
//
PolygonMesh square();

//          * \ .
//         /    e2
//      e3         \ .
//     /            *
//   *               /
//    \             e1
//     e4          /
//       \        *
//        *--e0-->
//
PolygonMesh pentagon();

//       -*<-e0-
//     /         \ .
//    |           |
//    |           |
//     \         /
//       -e1->*-
//
PolygonMesh digon();

//      --e0--
//    /        \ .
//   |          |
//   |          |
//    \        /
//      -->*--
//
PolygonMesh monogon();

//      -------
//    /         \ .
//   |  *-h0a->  |
//   |  <-h0b-*  |
//    \         /
//      -------
PolygonMesh bubble();

//      -------
//    /         \ .
//   |-->*--h0a--|
//   |--*<--h0b--|
//    \         /
//      -------
PolygonMesh dual_bubble();

//         5
//       /   \ .
//      / f1  \ .
//     /       \ .
//    2 -- 3 -- 4
//    |         |
//    |   f0    |
//    |         |
//    0 ------- 1
//
PolygonMesh glued_polygons();

//    6 ---- 7 ---- 8
//    |      |      |
//    |  f2  |  f3  |
//    3 ---- 4 ---- 5
//    |      |      |
//    |  f0  |  f1  |
//    0 ---- 1 ---- 2
//
PolygonMesh grid();

//    2 ----------- 3
//    |             |
//    |    6 - 7    |
//    |    |   |    |
//    |    4 - 5    |
//    |  /          |
//    0 ----------- 1
//
PolygonMesh annulus();

//    2 ---- 3   6 ---- 7
//    |      |   |      |
//    |  f0  |   |  f1  |
//    0 ---- 1   4 ---- 5
//
PolygonMesh two_squares();

//      <--h01--*
//    *           |
//    |           h10
//   h11          |
//    |           *
//      *--h01-->
//
PolygonMesh torus();

//      --e0--
//    /        \ .
//   |    ||    |
//   | h10||h11 |
//    \   ||   /
//      -->*--
//
PolygonMesh monogon_with_leaf();

// Generate a random polygon mesh with a given number of edges
PolygonMesh random_polygon_mesh(long num_edges, long rng_seed = 0);

} // namespace wmtk::tests
