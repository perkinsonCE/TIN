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
#include <cmath>

#include "globals.h"

// verifies if point coordinates are defined
bool point_not_undef(struct POINT p)
{
    if(p.x != UNDEF && p.y != UNDEF)
        return(true);
    else
        return(false);
}

// returns if 2 points have the same coordinates
bool points_are_equal(struct POINT p1, struct POINT p2)
{
    return(fabs(p1.x - p2.x) < epsVertex && fabs(p1.y - p2.y) < epsVertex);
}

// return if 2 edges have the same points
bool edges_are_equal(struct POINT p1, struct POINT p2, struct POINT p3, struct POINT p4)
{
    return((points_are_equal(p1, p3) == true && points_are_equal(p2, p4) == true) ||
           (points_are_equal(p1, p4) == true && points_are_equal(p2, p3) == true));
}

// verifies if points p1,p2 are an edge of a triangle adjacent to the parent triangle and returns the pointer to the adjacent triangle
//
// in addition, stores as pointers in a vector the coordinates of the point opposite to the edge and the points of the edge, in this order
struct TRIANGLE *adjacent_triangle_by_parent(struct TRIANGLE *parent, struct TRIANGLE *child, struct POINT p1, struct POINT p2, struct POINT *point_seq)
{
    struct TRIANGLE *t;
    int i, j;

    // process edges
    for(i = 0;i < 3;i++)
    {
        // adjacent triangle by edge
        t = parent->adjacent[i];

        // if adjacent triangle exists
        if(t != NULL)
        {
            // if points form an edge
            if(points_are_equal(p1, t->p_i) == true && points_are_equal(p2, t->p_j) == true)
            {
                // store points: opposite and from the edge
                point_seq[0] = t->p_k;
                point_seq[1] = t->p_i;
                point_seq[2] = t->p_j;
                
                // if necessary to change adjacency
                if(parent != child)
                    // define adjacent to adjacent triangle as child triangle
                    for(j = 0;j < 3;j++)
                        if(parent->adjacent[i]->adjacent[j] == parent)
                            parent->adjacent[i]->adjacent[j] = child;
                    
                return(t);
            }
            else
            // if points form an edge
            if(points_are_equal(p1, t->p_j) == true && points_are_equal(p2, t->p_i) == true)
            {
                // store points: opposite and from the edge
                point_seq[0] = t->p_k;
                point_seq[1] = t->p_j;
                point_seq[2] = t->p_i;
            
                // if necessary to change adjacency
                if(parent != child)
                    // define adjacent to adjacent triangle as child triangle
                    for(j = 0;j < 3;j++)
                        if(parent->adjacent[i]->adjacent[j] == parent)
                            parent->adjacent[i]->adjacent[j] = child;
                    
                return(t);
            }
            else
            // if points form an edge
            if(points_are_equal(p1, t->p_j) == true && points_are_equal(p2, t->p_k) == true)
            {
                // store points: opposite and from the edge
                point_seq[0] = t->p_i;
                point_seq[1] = t->p_j;
                point_seq[2] = t->p_k;
            
                // if necessary to change adjacency
                if(parent != child)
                    // define adjacent to adjacent triangle as child triangle
                    for(j = 0;j < 3;j++)
                        if(parent->adjacent[i]->adjacent[j] == parent)
                            parent->adjacent[i]->adjacent[j] = child;
                    
                return(t);
            }
            else
            // if points form an edge
            if(points_are_equal(p1, t->p_k) == true && points_are_equal(p2, t->p_j) == true)
            {
                // store points: opposite and from the edge
                point_seq[0] = t->p_i;
                point_seq[1] = t->p_k;
                point_seq[2] = t->p_j;
            
                // if necessary to change adjacency
                if(parent != child)
                    // define adjacent to adjacent triangle as child triangle
                    for(j = 0;j < 3;j++)
                        if(parent->adjacent[i]->adjacent[j] == parent)
                            parent->adjacent[i]->adjacent[j] = child;
                    
                return(t);
            }
            else
            // if points form an edge
            if(points_are_equal(p1, t->p_k) == true && points_are_equal(p2, t->p_i) == true)
            {
                // store points: opposite and from the edge
                point_seq[0] = t->p_j;
                point_seq[1] = t->p_k;
                point_seq[2] = t->p_i;
            
                // if necessary to change adjacency
                if(parent != child)
                    // define adjacent to adjacent triangle as child triangle
                    for(j = 0;j < 3;j++)
                        if(parent->adjacent[i]->adjacent[j] == parent)
                            parent->adjacent[i]->adjacent[j] = child;
                    
                return(t);
            }
            else
            // if points form an edge
            if(points_are_equal(p1, t->p_i) == true && points_are_equal(p2, t->p_k) == true)
            {
                // store points: opposite and from the edge
                point_seq[0] = t->p_j;
                point_seq[1] = t->p_i;
                point_seq[2] = t->p_k;
            
                // if necessary to change adjacency
                if(parent != child)
                    // define adjacent to adjacent triangle as child triangle
                    for(j = 0;j < 3;j++)
                        if(parent->adjacent[i]->adjacent[j] == parent)
                            parent->adjacent[i]->adjacent[j] = child;
                    
                return(t);
            }
        }
    }
    
    // no adjacent triangle by the given edge
    t = NULL;
    
    return(t);
}

// verifies if points p1,p2 form an edge of triangle t and return the opposite point
struct POINT opposite_point_by_edge(struct POINT p1, struct POINT p2, struct TRIANGLE *t)
{
    struct POINT opposite_point;

    // if triangle exists
    if(t != NULL)
    {
        // verifies edge i,j
        if(edges_are_equal(p1, p2, t->p_i, t->p_j) == true)
        {
            opposite_point = t->p_k;
                
            return(opposite_point);
        }
        // verifies edge j,k
        if(edges_are_equal(p1, p2, t->p_j, t->p_k) == true)
        {
            opposite_point = t->p_i;
                
            return(opposite_point);
        }
        // verifies edge i,k
        if(edges_are_equal(p1, p2, t->p_i, t->p_k) == true)
        {
            opposite_point = t->p_j;
                
            return(opposite_point);
        }
    }
    
    // set opposite point with undefined coordinates as points p1,p2 are not an edge of triangle t
    opposite_point.x = UNDEF;
    opposite_point.y = UNDEF;
    opposite_point.z = UNDEF;
    opposite_point.address = NULL;

    return(opposite_point);
}

