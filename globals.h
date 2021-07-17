/*
    Copyright(c) 2021

    This file is part of an implementation of the Delaunay triangulation
	written in C++ for standard C++ compilers. Please, feel free to
	redistribute and/or modify it, but without any warrant.
	
	The code was written as part of the Master's dissertation by
	Henrique Rennó de Azeredo Freitas <henrique.renno@inpe.br>, which is
	available from http://urlib.net/8JMKD3MGP3W/3HDPRRS. Any doubt concerning
	the code can be asked to this author.
	
	If you use this code to publish other related works, please cite
	the appropriate reference(s):
	
	Freitas, H. R. A., 2014. Drainage networks and watersheds delineation
	derived from TIN-based digital elevation models. Master's dissertation.
	Instituto Nacional de Pesquisas Espaciais (INPE). Available from:
	http://urlib.net/8JMKD3MGP3W/3HDPRRS
	
	Freitas, H. R. A., 2015. Description of a C++ implementation of the
	incremental algorithm to generate the Delaunay triangulation.
	Available from:	http://urlib.net/8JMKD3MGP3W34P/3JJUR2L
*/

#include <vector>

#include "structs.h"

using namespace std;

// points vector
extern vector <struct POINT> points;
// 3 points of enclosing triangle
extern struct POINT p_0, p_1, p_2;
// root node of triangles tree
extern struct TRIANGLE *triangles_tree;
// flags to indicate if point is on an edge
extern bool same_slope_ij, same_slope_jk, same_slope_ki;
// Delaunay triangles vector
extern vector <struct TRIANGLE *> delaunay_triangles;
// error value for vertices
const double epsVertex = 1.e-5;
// error value for edges
const double epsEdge = 1.e-6;
// constant undefined value
const double UNDEF = -9999.;
// increment in outermost coordinates
const double INC = 10.;
