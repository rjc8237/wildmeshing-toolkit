#pragma once

#include <wmtk/PolygonMesh.hpp>

namespace wmtk::tests {


//        * \ 
//       /   e1
//      e2    \ 
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

//			    * \ 
//			   /  	e2
//		  e3   	  	\ 
//		 /			      *
//   *               /
//    \             e1
//     e4          /
//       \        *
//        *--e0-->
//
PolygonMesh pentagon();

//       -*<-e0-
//     /         \ 
//    |           |
//    |           |
//     \         /
//       -e1->*-
//
PolygonMesh digon();

//      --e0--
//    /        \ 
//   |          |
//   |          |
//    \        /
//      -->*--
//
PolygonMesh monogon();

//      -------
//    /         \ 
//   |  *-h0a->  |
//   |  <-h0b-*  |
//    \         /
//      -------
PolygonMesh bubble();

//      -------
//    /         \ 
//   |-->*--h0a--|
//   |--*<--h0b--|
//    \         /
//      -------
PolygonMesh dual_bubble();

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
PolygonMesh grid();


// TODO We could also use any TriMesh example as input. This could be made explicit in
// this file, but I am hesitant to duplicate the ascii art in case it is changed later.

} // namespace wmtk::tests
