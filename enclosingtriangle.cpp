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

#include <cstdio>
#include <cstdlib>

#include "globals.h"

// 3 points of enclosing triangle
struct POINT p_0, p_1, p_2;
// root node of triangles tree
struct TRIANGLE *triangles_tree;

// define the coordinates of the points p0,p1,p2 of the
// enclosing triangle of a set of points p in the plane
//
// the coordinates are calculated relative to an
// enclosing rectangle
//
// with the maximum point (max_x+INC,max_y+INC) of the
// rectangle defined, the points (0,K) and (K,0) of the
// triangle are calculated from the intersection of the
// line with slope = -1 (-45o) that passes through the
// maximum point and intersects the x and y axes.
// the last point is defined with coordinates (-K,-K)
void define_enclosing_triangle(vector <struct POINT> &p)
{
    int i;
    double max_x, max_y, K, INC;
    struct TRIANGLE t;

    // define the increment in coordinates
    INC = 10.;
        
    // define maximum coordinates in x and y
    max_x = -1.;
    max_y = -1.;
    for(i = 0;i < p.size();i++)
    {
        if(p[i].x > max_x)
            max_x = p[i].x;
        if(p[i].y > max_y)
            max_y = p[i].y;
    }
    
    // define maximum coordinates of the enclosing
    // rectangle incremented by INC shifting up and right
    max_x += INC;
    max_y += INC;
    
    // define coordinates of points p0,p1,p2 with
    // K = max_x + max_y given by the equation
    // (K - max_y)/(0 - max_x) = -1 of the line that passes
    // through the points (0,K) and (max_x,max_y) with
    // slope = -1
    K = max_x + max_y;
    // multiply K by a constant
    K *= INC;
  
    // coordinates of points p0,p1,p2
    p_0.x = -K;
    p_0.y = -K;
    p_1.x = 0.;
    p_1.y = K;
    p_2.x = K;
    p_2.y = 0.;
    
    // define references of the enclosing triangle points
    p_0.address = &p_0;
    p_1.address = &p_1;
    p_2.address = &p_2;
    
    // define enclosing triangle p0,p1,p2
    t.p_i = p_0;
    t.p_j = p_1;
    t.p_k = p_2;
    t.number_of_children = 0;
    t.children = NULL;
    t.number_of_parents = 0;
    for(i = 0;i < 2;i++)
        t.parents[i] = NULL;
    for(i = 0;i < 3;i++)
        t.adjacent[i] = NULL;
    t.triangle_index = -1;
    t.processed = false;
    t.visited = false;
    
    // allocate memory to reference a tree of triangles
    triangles_tree = (struct TRIANGLE *) malloc(sizeof(struct TRIANGLE));
    // define triangle p0,p1,p2 as child of first tree node
    triangles_tree->children = (struct TRIANGLE **) malloc(sizeof(struct TRIANGLE *));
    triangles_tree->children[0] = (struct TRIANGLE *) malloc(sizeof(struct TRIANGLE));
    *(triangles_tree->children[0]) = t;
    triangles_tree->number_of_children = 1;
    triangles_tree->processed = false;
    triangles_tree->visited = false;

    return;
}
