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
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <queue>

#include "globals.h"
#include "functions.h"

// returns flag indicating if triangle has p0, p1 or p2 as vertex
bool valid_delaunay_triangle(struct TRIANGLE *t)
{
    if(t == NULL)
        return(false);
        
    if((points_are_equal(t->p_i, p_0) == true) ||
       (points_are_equal(t->p_i, p_1) == true) ||
       (points_are_equal(t->p_i, p_2) == true) ||
       (points_are_equal(t->p_j, p_0) == true) ||
       (points_are_equal(t->p_j, p_1) == true) ||
       (points_are_equal(t->p_j, p_2) == true) ||
       (points_are_equal(t->p_k, p_0) == true) ||
       (points_are_equal(t->p_k, p_1) == true) ||
       (points_are_equal(t->p_k, p_2) == true))
        return(false);
    
    return(true);
}


// perform a depth-first search to select the Delaunay triangles from the leaf nodes of the triangles tree
void define_delaunay_triangles_from_tree(struct TRIANGLE *t, vector <struct TRIANGLE *> &triangles)
{
    int i;

    // if triangle is not visited
    if(t->visited == false)
    {
        // sets triangle as visited
        t->visited = true;

        // if triangle has children is not leaf, search continues
        if(t->number_of_children > 0)
        {
            for(i = 0;i < t->number_of_children;i++)
                define_delaunay_triangles_from_tree(t->children[i], triangles);
        }
        // if triangle has no children, then it is a Delaunay triangle
        else
        {
            // verifies if triangle does not include points p0,p1,p2 from enclosing triangle
            if(valid_delaunay_triangle(t) == true)
            {
                // inserts triangle in the Delaunay triangles vector
                triangles.push_back(t);
	    }
        }
    }
    
    return;
}

// an edge p1,p2 is legalized by calculating the center and radius of the circunference from points p,p1,p2 and by verifing if the point of the adjacent triangle opposite to the edge is inside the circunference, so the edge is illegal and it is rotated creating 2 new triangles
void legalize_edge(struct POINT p, struct POINT p1, struct POINT p2, struct TRIANGLE *t_o, struct TRIANGLE *t_p)
{
    double xc, yc, r, dist;
    struct POINT opposite_point, center, useless_point_seq[3];
    int i;

    // cases where some points should be excluded from the set
    //
    // if point p is equal to edge point
    if(points_are_equal(p, p1) == true)
    {
        cout << "point p equals to p1" << endl;
        cout << "p.index = " << p.point_index << endl;
        cout << "p1.index = " << p1.point_index << endl;
        cout << setprecision(15) << "p.x = " << p.x << " | p.y = " << p.y << " | p.z = " << p.z << endl;
        cout << setprecision(15) << "p1.x = " << p1.x << " | p1.y = " << p1.y << " | p1.z = " << p1.z << endl;
        exit(1);
        return;
    }
    // if point p is equal to edge point
    if(points_are_equal(p, p2) == true)
    {
        cout << "point p equals to p2" << endl;
        cout << "p.index = " << p.point_index << endl;
        cout << "p2.index = " << p2.point_index << endl;
        cout << setprecision(15) << "p.x = " << p.x << " | p.y = " << p.y << " | p.z = " << p.z << endl;
        cout << setprecision(15) << "p2.x = " << p2.x << " | p2.y = " << p2.y << " | p2.z = " << p2.z << endl;
        exit(1);
        return;
    }
    // if points p1 and p2 are equal
    if(points_are_equal(p1, p2) == true)
    {
        cout << "point p1 equals to p2" << endl;
        cout << "p1.index = " << p1.point_index << endl;
        cout << "p2.index = " << p2.point_index << endl;
        cout << setprecision(15) << "p1.x = " << p1.x << " | p1.y = " << p1.y << " | p1.z = " << p1.z << endl;
        cout << setprecision(15) << "p2.x = " << p2.x << " | p2.y = " << p2.y << " | p2.z = " << p2.z << endl;
        exit(1);
        return;
    }
    // if p,p1,p2 are colinear in x the circunference does not exist
    if(fabs(p.x - p1.x) < epsVertex && fabs(p1.x - p2.x) < epsVertex)
    {
        cout << "3 colinear points, vertical line" << endl;
        cout << setprecision(15) << "p.x =  " << p.point_index << " " << p.x <<  " | p.y = " << p.y << endl;
        cout << setprecision(15) << "p1.x = " << p1.point_index << " " << p1.x << " | p1.y = " << p1.y << endl;
        cout << setprecision(15) << "p2.x = " << p2.point_index << " " << p2.x << " | p2.y = " << p2.y << endl;
        exit(1);
    }
    // if p,p1,p2 are colinear in y the circunference does not exist
    if(fabs(p.y - p1.y) < epsVertex && fabs(p1.y - p2.y) < epsVertex)
    {
        cout << "3 colinear points, horizontal line" << endl;
        cout << setprecision(15) << "p.x =  " << p.point_index << " " << p.x <<  " | p.y = " << p.y << endl;
        cout << setprecision(15) << "p1.x = " << p1.point_index << " " << p1.x << " | p1.y = " << p1.y << endl;
        cout << setprecision(15) << "p2.x = " << p2.point_index << " " << p2.x << " | p2.y = " << p2.y << endl;
        exit(1);
    }
    // if p,p1,p2 are colinear with same slope the circunference does not exist
    if(same_slope(p, p1, p2) == true)
    {
        cout << "3 colinear points, same slope" << endl;
        cout << setprecision(15) << "p.x =  " << p.x <<  " | p.y = " << p.y << endl;
        cout << setprecision(15) << "p1.x = " << p1.x << " | p1.y = " << p1.y << endl;
        cout << setprecision(15) << "p2.x = " << p2.x << " | p2.y = " << p2.y << endl;
        exit(1);
    }
    
    // define point opposite to the edge p1,p2 from adjacent triangle
    opposite_point = opposite_point_by_edge(p1, p2, t_o);
    // if point is not defined
    if(point_not_undef(opposite_point) == false)
        return;

    // calculates center and radius of circunference passing through points p,p1,p2
    circle_by_3_points(p, p1, p2, &xc, &yc, &r);

    // verifies if point opposite to the edge p1,p2 is inside the circunference
    //
    // calculates distance from circunference center to the opposite point
    center.x = xc;
    center.y = yc;
    dist = distance_between(center, opposite_point);

    // verifies if distance is less than radius, so the edge is not valid
    if(dist < r)
    {
            // divide current triangles adding new triangles p,p1,opposite and p,p2,opposite as children
            //
            // triangle with point p
            t_p->children = (struct TRIANGLE **) malloc(2*sizeof(struct TRIANGLE *));
            t_p->number_of_children = 2;
            for(i = 0;i < t_p->number_of_children;i++)
                t_p->children[i] = (struct TRIANGLE *) malloc(sizeof(struct TRIANGLE));
            
            // triangle p,p1,opposite
            t_p->children[0]->p_i = p;
            t_p->children[0]->p_j = p1;
            t_p->children[0]->p_k = opposite_point;
            t_p->children[0]->children = NULL;
            t_p->children[0]->number_of_children = 0;
            t_p->children[0]->parents[0] = t_p;
            t_p->children[0]->parents[1] = t_o;
            t_p->children[0]->number_of_parents = 2;
            t_p->children[0]->processed = false;
            t_p->children[0]->visited = false;
            // triangle p,p2,opposite
            t_p->children[1]->p_i = p;
            t_p->children[1]->p_j = p2;
            t_p->children[1]->p_k = opposite_point;
            t_p->children[1]->children = NULL;
            t_p->children[1]->number_of_children = 0;
            t_p->children[1]->parents[0] = t_p;
            t_p->children[1]->parents[1] = t_o;
            t_p->children[1]->number_of_parents = 2;
            t_p->children[1]->processed = false;
            t_p->children[1]->visited = false;

            // define adjacency between new triangles
            //
            // triangle p,p1,opposite
            // adjacent by parent p1,p2,p
            t_p->children[0]->adjacent[0] = adjacent_triangle_by_parent(t_p, t_p->children[0], p1, p, useless_point_seq);
            // adjacent by parent p1,p2,opposite
            t_p->children[0]->adjacent[1] = adjacent_triangle_by_parent(t_o, t_p->children[0], p1, opposite_point, useless_point_seq);
            // adjacent p,p2,opposite
            t_p->children[0]->adjacent[2] = t_p->children[1];
            //
            // triangle p,p2,opposite
            // adjacente pelo triângulo pai p1,p2,p
            t_p->children[1]->adjacent[0] = adjacent_triangle_by_parent(t_p, t_p->children[1], p2, p, useless_point_seq);
            // adjacent by parent p1,p2,opposite
            t_p->children[1]->adjacent[1] = adjacent_triangle_by_parent(t_o, t_p->children[1], p2, opposite_point, useless_point_seq);
            // adjacent p,p1,opposite
            t_p->children[1]->adjacent[2] = t_p->children[0];

            // triangle with opposite point, children are the same of triangle with point p
            t_o->children = (struct TRIANGLE **) malloc(2*sizeof(struct TRIANGLE *));
            t_o->number_of_children = 2;

            // triangle p,p1,opposite
            t_o->children[0] = t_p->children[0];
            // triangle p,p2,opposite
            t_o->children[1] = t_p->children[1];
            
            // recursively legalize edges p1,opposite e opposite,p2
            //
            // adjacent triangles of new triangles
            // adjacent by edge p1,opposite
            legalize_edge(p, p1, opposite_point, t_o->children[0]->adjacent[1], t_o->children[0]);
            // adjacent by edge p2,opposite
            legalize_edge(p, opposite_point, p2, t_o->children[1]->adjacent[1], t_o->children[1]);
    }

    return;
}

// calculates the Delaunay triangulation from a set of points with the incremental algorithm
//
// the points are inserted one at a time and the triangulation is locally modified with the insertion of each point
vector <struct TRIANGLE *> delaunay_triangulation(vector <struct POINT> &p)
{
    vector <struct TRIANGLE *> triangles;
    queue <struct TRIANGLE *> tq;
    struct TRIANGLE *t_node, *t_adjacent;
    struct POINT point_seq[3], useless_point_seq[3];
    int point_in_triangle_result, p_index, t_index;
    int i;
    
    // inserts a point and verifies the triangle that contains the point
    for(p_index = 0;p_index < p.size();p_index++)
    {
        // search the triangle that contains the point traversing the tree
        //
        // sets the enclosing triangle as not processed and inserts it into the queue
        triangles_tree->children[0]->processed = false;
        tq.push(triangles_tree->children[0]);

        // while the queue is not empty, triangles are still verified
        while(tq.empty() == false)
        {
            // remove triangle from queue
            t_node = tq.front();
            tq.pop();
            
            // verifies if triangle is already processed, which can occur when the point is located exactly on an edge
            //
            // if triangle is already processed, continue, otherwise sets it as processed
            if(t_node->processed == true)
                continue;
            else
                t_node->processed = true;
            
            // get location of point in the triangle
            point_in_triangle_result = point_in_triangle(p[p_index], *t_node);

            // if the triangle has children, add the children into the queue if the point is on the triangle
            if(t_node->number_of_children > 0)
            {
                // if point is inside or on an edge of the triangle
                if(point_in_triangle_result > 0)
                {
                    // no need to verify other triangles in the queue, empty the queue
                    tq = queue<struct TRIANGLE *>();

                    // add children in the queue
                    for(t_index = 0;t_index < t_node->number_of_children;t_index++)
                    {
                        // sets triangle as not processed
                        t_node->children[t_index]->processed = false;
                        tq.push(t_node->children[t_index]);
                    }
                }
            }
            // if triangle is leaf node, divide triangle according to location of point in triangle
            else
            {
                // if point is inside triangle
                if(point_in_triangle_result == 2)
                {
                    // add edges from the point to the vertices of the triangle dividing it into 3 new triangles
                    t_node->children = (struct TRIANGLE **) malloc(3*sizeof(struct TRIANGLE *));
                    t_node->number_of_children = 3;
                    for(i = 0;i < t_node->number_of_children;i++)
                        t_node->children[i] = (struct TRIANGLE *) malloc(sizeof(struct TRIANGLE));
                    
                    // triangle i,j,p
                    t_node->children[0]->p_i = t_node->p_i;
                    t_node->children[0]->p_j = t_node->p_j;
                    t_node->children[0]->p_k = p[p_index];
                    t_node->children[0]->children = NULL;
                    t_node->children[0]->number_of_children = 0;
                    t_node->children[0]->parents[0] = t_node;
                    t_node->children[0]->number_of_parents = 1;
                    t_node->children[0]->processed = false;
                    t_node->children[0]->visited = false;
                    // triangle j,k,p
                    t_node->children[1]->p_i = t_node->p_j;
                    t_node->children[1]->p_j = t_node->p_k;
                    t_node->children[1]->p_k = p[p_index];
                    t_node->children[1]->children = NULL;
                    t_node->children[1]->number_of_children = 0;
                    t_node->children[1]->parents[0] = t_node;
                    t_node->children[1]->number_of_parents = 1;
                    t_node->children[1]->processed = false;
                    t_node->children[1]->visited = false;
                    // triangle i,k,p
                    t_node->children[2]->p_i = t_node->p_i;
                    t_node->children[2]->p_j = t_node->p_k;
                    t_node->children[2]->p_k = p[p_index];
                    t_node->children[2]->children = NULL;
                    t_node->children[2]->number_of_children = 0;
                    t_node->children[2]->parents[0] = t_node;
                    t_node->children[2]->number_of_parents = 1;
                    t_node->children[2]->processed = false;
                    t_node->children[2]->visited = false;
                    
                    // define adjacency between new triangles
                    //
                    // triangle i,j,p:
                    // adjacent by parent
                    t_node->children[0]->adjacent[0] = adjacent_triangle_by_parent(t_node, t_node->children[0], t_node->p_i, t_node->p_j, useless_point_seq);
                    // adjacent j,k,p
                    t_node->children[0]->adjacent[1] = t_node->children[1];
                    // adjacent i,k,p
                    t_node->children[0]->adjacent[2] = t_node->children[2];
                    //
                    // triangle j,k,p:
                    // adjacent by parent
                    t_node->children[1]->adjacent[0] = adjacent_triangle_by_parent(t_node, t_node->children[1], t_node->p_j, t_node->p_k, useless_point_seq);
                    // adjacent i,k,p
                    t_node->children[1]->adjacent[1] = t_node->children[2];
                    // adjacent i,j,p
                    t_node->children[1]->adjacent[2] = t_node->children[0];
                    //
                    // triangle i,k,p:
                    // adjacent by parent
                    t_node->children[2]->adjacent[0] = adjacent_triangle_by_parent(t_node, t_node->children[2], t_node->p_i, t_node->p_k, useless_point_seq);
                    // adjacent j,k,p
                    t_node->children[2]->adjacent[1] = t_node->children[1];
                    // adjacent i,j,p
                    t_node->children[2]->adjacent[2] = t_node->children[0];

                    // legalize the edges of the divided triangle
                    legalize_edge(p[p_index], t_node->p_i, t_node->p_j, t_node->children[0]->adjacent[0], t_node->children[0]);
                    legalize_edge(p[p_index], t_node->p_j, t_node->p_k, t_node->children[1]->adjacent[0], t_node->children[1]);
                    legalize_edge(p[p_index], t_node->p_k, t_node->p_i, t_node->children[2]->adjacent[0], t_node->children[2]);
                    
                    // no need to verify other triangles in the queue, empty the queue
                    tq = queue<struct TRIANGLE *>();
                }
                // if the point is on an edge
                else
                if(point_in_triangle_result == 1)
                {
                    // add an edge from the point to the point of the triangle opposite to the edge dividing it into 2 new triangles
                    t_node->children = (struct TRIANGLE **) malloc(2*sizeof(struct TRIANGLE *));
                    t_node->number_of_children = 2;
                    for(i = 0;i < t_node->number_of_children;i++)
                        t_node->children[i] = (struct TRIANGLE *) malloc(sizeof(struct TRIANGLE));
                    
                    // verifies on what edge the point is located and defines the new triangles with the opposite point
                    //
                    // edge i,j
                    // adjacent triangle: point l <-> point k
                    if(same_slope_ij == true)
                    {
                        // triangle i,k,p
                        t_node->children[0]->p_i = t_node->p_i;
                        t_node->children[0]->p_j = t_node->p_k;
                        t_node->children[0]->p_k = p[p_index];
                        t_node->children[0]->children = NULL;
                        t_node->children[0]->number_of_children = 0;
                        t_node->children[0]->parents[0] = t_node;
                        t_node->children[0]->number_of_parents = 1;
                        t_node->children[0]->processed = false;
                        t_node->children[0]->visited = false;
                        // triangle j,k,p
                        t_node->children[1]->p_i = t_node->p_j;
                        t_node->children[1]->p_j = t_node->p_k;
                        t_node->children[1]->p_k = p[p_index];
                        t_node->children[1]->children = NULL;
                        t_node->children[1]->number_of_children = 0;                            
                        t_node->children[1]->parents[0] = t_node;
                        t_node->children[1]->number_of_parents = 1;
                        t_node->children[1]->processed = false;
                        t_node->children[1]->visited = false;

                        // define 2 new triangles from the adjacent triangle by the edge
                        t_adjacent = adjacent_triangle_by_parent(t_node, t_node, t_node->p_i, t_node->p_j, point_seq);
                        
                        if(t_adjacent != NULL)
                        {
                            t_adjacent->children = (struct TRIANGLE **) malloc(2*sizeof(struct TRIANGLE *));
                            t_adjacent->number_of_children = 2;
                            for(i = 0;i < t_adjacent->number_of_children;i++)
                                t_adjacent->children[i] = (struct TRIANGLE *) malloc(sizeof(struct TRIANGLE));

                            // sets adjacent triangle as processed
                            t_adjacent->processed = true;
                            
                            // triangle i,l,p
                            t_adjacent->children[0]->p_i = point_seq[1];
                            t_adjacent->children[0]->p_j = point_seq[0];
                            t_adjacent->children[0]->p_k = p[p_index];
                            t_adjacent->children[0]->children = NULL;
                            t_adjacent->children[0]->number_of_children = 0;
                            t_adjacent->children[0]->parents[0] = t_adjacent;
                            t_adjacent->children[0]->number_of_parents = 1;
                            t_adjacent->children[0]->processed = false;
                            t_adjacent->children[0]->visited = false;
                            // triangle j,l,p
                            t_adjacent->children[1]->p_i = point_seq[2];
                            t_adjacent->children[1]->p_j = point_seq[0];
                            t_adjacent->children[1]->p_k = p[p_index];
                            t_adjacent->children[1]->children = NULL;
                            t_adjacent->children[1]->number_of_children = 0;
                            t_adjacent->children[1]->parents[0] = t_adjacent;
                            t_adjacent->children[1]->number_of_parents = 1;
                            t_adjacent->children[1]->processed = false;
                            t_adjacent->children[1]->visited = false;
                        }
                        
                        // define adjacency between new triangles
                        //
                        // triangle i,k,p:
                        // adjacent by parent
                        t_node->children[0]->adjacent[0] = adjacent_triangle_by_parent(t_node, t_node->children[0], t_node->p_i, t_node->p_k, useless_point_seq);
                        // adjacent j,k,p
                        t_node->children[0]->adjacent[1] = t_node->children[1];
                        // adjacent i,l,p
                        if(t_adjacent != NULL)
                            t_node->children[0]->adjacent[2] = t_adjacent->children[0];

                        // triangle j,k,p:
                        // adjacent by parent
                        t_node->children[1]->adjacent[0] = adjacent_triangle_by_parent(t_node, t_node->children[1], t_node->p_j, t_node->p_k, useless_point_seq);
                        // adjacent i,k,p
                        t_node->children[1]->adjacent[1] = t_node->children[0];
                        // adjacent j,l,p
                        if(t_adjacent != NULL)
                            t_node->children[1]->adjacent[2] = t_adjacent->children[1];
                        
                        if(t_adjacent != NULL)
                        {
                            // triangle i,l,p:
                            // adjacente pelo triângulo pai
                            t_adjacent->children[0]->adjacent[0] = adjacent_triangle_by_parent(t_adjacent, t_adjacent->children[0], point_seq[1], point_seq[0], useless_point_seq);
                            // adjacent j,l,p
                            t_adjacent->children[0]->adjacent[1] = t_adjacent->children[1];
                            // adjacent i,k,p
                            t_adjacent->children[0]->adjacent[2] = t_node->children[0];
                            //
                            // triangle j,l,p:
                            // adjacente pelo triângulo pai
                            t_adjacent->children[1]->adjacent[0] = adjacent_triangle_by_parent(t_adjacent, t_adjacent->children[1], point_seq[2], point_seq[0], useless_point_seq);
                            // adjacent i,l,p
                            t_adjacent->children[1]->adjacent[1] = t_adjacent->children[0];
                            // adjacent j,k,p
                            t_adjacent->children[1]->adjacent[2] = t_node->children[1];
                        }

                        // legalize the edges of the divided triangle
                        legalize_edge(p[p_index], t_node->p_i, point_seq[0], t_adjacent->children[0]->adjacent[0], t_adjacent->children[0]);
                        legalize_edge(p[p_index], point_seq[0], t_node->p_j, t_adjacent->children[1]->adjacent[0], t_adjacent->children[1]);
                        legalize_edge(p[p_index], t_node->p_j, t_node->p_k, t_node->children[1]->adjacent[0], t_node->children[1]);
                        legalize_edge(p[p_index], t_node->p_k, t_node->p_i, t_node->children[0]->adjacent[0], t_node->children[0]);
                    }
                    else
                    // edge j,k
                    // adjacent triangle: point l <-> point i
                    if(same_slope_jk == true)
                    {
                        // triangle j,i,p
                        t_node->children[0]->p_i = t_node->p_j;
                        t_node->children[0]->p_j = t_node->p_i;
                        t_node->children[0]->p_k = p[p_index];
                        t_node->children[0]->children = NULL;
                        t_node->children[0]->number_of_children = 0;
                        t_node->children[0]->parents[0] = t_node;
                        t_node->children[0]->number_of_parents = 1;
                        t_node->children[0]->processed = false;
                        t_node->children[0]->visited = false;
                        // triangle k,i,p
                        t_node->children[1]->p_i = t_node->p_k;
                        t_node->children[1]->p_j = t_node->p_i;
                        t_node->children[1]->p_k = p[p_index];
                        t_node->children[1]->children = NULL;
                        t_node->children[1]->number_of_children = 0;
                        t_node->children[1]->parents[0] = t_node;
                        t_node->children[1]->number_of_parents = 1;
                        t_node->children[1]->processed = false;
                        t_node->children[1]->visited = false;
                        
                        // define 2 new triangles from the adjacent triangle by the edge
                        t_adjacent = adjacent_triangle_by_parent(t_node, t_node, t_node->p_j, t_node->p_k, point_seq);

                        if(t_adjacent != NULL)
                        {
                            t_adjacent->children = (struct TRIANGLE **) malloc(2*sizeof(struct TRIANGLE *));
                            t_adjacent->number_of_children = 2;
                            for(i = 0;i < t_adjacent->number_of_children;i++)
                                t_adjacent->children[i] = (struct TRIANGLE *) malloc(sizeof(struct TRIANGLE));

                            // sets adjacent triangle as processed
                            t_adjacent->processed = true;
                            
                            // triangle j,l,p
                            t_adjacent->children[0]->p_i = point_seq[1];
                            t_adjacent->children[0]->p_j = point_seq[0];
                            t_adjacent->children[0]->p_k = p[p_index];
                            t_adjacent->children[0]->children = NULL;
                            t_adjacent->children[0]->number_of_children = 0;
                            t_adjacent->children[0]->parents[0] = t_adjacent;
                            t_adjacent->children[0]->number_of_parents = 1;
                            t_adjacent->children[0]->processed = false;
                            t_adjacent->children[0]->visited = false;
                            // triangle k,l,p
                            t_adjacent->children[1]->p_i = point_seq[2];
                            t_adjacent->children[1]->p_j = point_seq[0];
                            t_adjacent->children[1]->p_k = p[p_index];
                            t_adjacent->children[1]->children = NULL;
                            t_adjacent->children[1]->number_of_children = 0;
                            t_adjacent->children[1]->parents[0] = t_adjacent;
                            t_adjacent->children[1]->number_of_parents = 1;
                            t_adjacent->children[1]->processed = false;
                            t_adjacent->children[1]->visited = false;
                        }
                        
                        // define adjacency between new triangles
                        //
                        // triangle j,i,p:
                        // adjacent by parent
                        t_node->children[0]->adjacent[0] = adjacent_triangle_by_parent(t_node, t_node->children[0], t_node->p_j, t_node->p_i, useless_point_seq);
                        // adjacent k,i,p
                        t_node->children[0]->adjacent[1] = t_node->children[1];
                        // adjacent j,l,p
                        if(t_adjacent != NULL)
                            t_node->children[0]->adjacent[2] = t_adjacent->children[0];

                        // triângulo k,i,p:
                        // adjacent by parent
                        t_node->children[1]->adjacent[0] = adjacent_triangle_by_parent(t_node, t_node->children[1], t_node->p_k, t_node->p_i, useless_point_seq);
                        // adjacent j,i,p
                        t_node->children[1]->adjacent[1] = t_node->children[0];
                        // adjacent k,l,p
                        if(t_adjacent != NULL)
                            t_node->children[1]->adjacent[2] = t_adjacent->children[1];
                        
                        if(t_adjacent != NULL)
                        {
                            // triangle j,l,p:
                            // adjacent by parent
                            t_adjacent->children[0]->adjacent[0] = adjacent_triangle_by_parent(t_adjacent, t_adjacent->children[0], point_seq[1], point_seq[0], useless_point_seq);
                            // adjacent k,l,p
                            t_adjacent->children[0]->adjacent[1] = t_adjacent->children[1];
                            // adjacent j,i,p
                            t_adjacent->children[0]->adjacent[2] = t_node->children[0];
                            //
                            // triangle k,l,p:
                            // adjacente pelo triângulo pai
                            t_adjacent->children[1]->adjacent[0] = adjacent_triangle_by_parent(t_adjacent, t_adjacent->children[1], point_seq[2], point_seq[0], useless_point_seq);
                            // adjacent j,l,p
                            t_adjacent->children[1]->adjacent[1] = t_adjacent->children[0];
                            // adjacent k,i,p
                            t_adjacent->children[1]->adjacent[2] = t_node->children[1];
                        }

                        // legalize the edges of the divided triangle
                        legalize_edge(p[p_index], t_node->p_j, point_seq[0], t_adjacent->children[0]->adjacent[0], t_adjacent->children[0]);
                        legalize_edge(p[p_index], point_seq[0], t_node->p_k, t_adjacent->children[1]->adjacent[0], t_adjacent->children[1]);
                        legalize_edge(p[p_index], t_node->p_k, t_node->p_i, t_node->children[1]->adjacent[0], t_node->children[1]);
                        legalize_edge(p[p_index], t_node->p_i, t_node->p_j, t_node->children[0]->adjacent[0], t_node->children[0]);
                    }
                    else
                    // edge k,i
                    // adjacent triangle: point l <-> point j
                    if(same_slope_ki == true)
                    {
                        // triangle k,j,p
                        t_node->children[0]->p_i = t_node->p_k;
                        t_node->children[0]->p_j = t_node->p_j;
                        t_node->children[0]->p_k = p[p_index];
                        t_node->children[0]->children = NULL;
                        t_node->children[0]->number_of_children = 0;
                        t_node->children[0]->parents[0] = t_node;
                        t_node->children[0]->number_of_parents = 1;
                        t_node->children[0]->processed = false;
                        t_node->children[0]->visited = false;
                        // triangle i,j,p
                        t_node->children[1]->p_i = t_node->p_i;
                        t_node->children[1]->p_j = t_node->p_j;
                        t_node->children[1]->p_k = p[p_index];
                        t_node->children[1]->children = NULL;
                        t_node->children[1]->number_of_children = 0;
                        t_node->children[1]->parents[0] = t_node;
                        t_node->children[1]->number_of_parents = 1;
                        t_node->children[1]->processed = false;
                        t_node->children[1]->visited = false;

                        // define 2 new triangles from the adjacent triangle by the edge
                        t_adjacent = adjacent_triangle_by_parent(t_node, t_node, t_node->p_k, t_node->p_i, point_seq);

                        if(t_adjacent != NULL)
                        {
                            t_adjacent->children = (struct TRIANGLE **) malloc(2*sizeof(struct TRIANGLE *));
                            t_adjacent->number_of_children = 2;
                            for(i = 0;i < t_adjacent->number_of_children;i++)
                                t_adjacent->children[i] = (struct TRIANGLE *) malloc(sizeof(struct TRIANGLE));

                            // sets adjacent triangle as processed
                            t_adjacent->processed = true;
                            
                            // triangle k,l,p
                            t_adjacent->children[0]->p_i = point_seq[1];
                            t_adjacent->children[0]->p_j = point_seq[0];
                            t_adjacent->children[0]->p_k = p[p_index];
                            t_adjacent->children[0]->children = NULL;
                            t_adjacent->children[0]->number_of_children = 0;
                            t_adjacent->children[0]->parents[0] = t_adjacent;
                            t_adjacent->children[0]->number_of_parents = 1;
                            t_adjacent->children[0]->processed = false;
                            t_adjacent->children[0]->visited = false;
                            // triangle i,l,p
                            t_adjacent->children[1]->p_i = point_seq[2];
                            t_adjacent->children[1]->p_j = point_seq[0];
                            t_adjacent->children[1]->p_k = p[p_index];
                            t_adjacent->children[1]->children = NULL;
                            t_adjacent->children[1]->number_of_children = 0;
                            t_adjacent->children[1]->parents[0] = t_adjacent;
                            t_adjacent->children[1]->number_of_parents = 1;
                            t_adjacent->children[1]->processed = false;
                            t_adjacent->children[1]->visited = false;
                        }
                        
                        // define adjacency between new triangles
                        //
                        // triangle k,j,p:
                        // adjacent by parent
                        t_node->children[0]->adjacent[0] = adjacent_triangle_by_parent(t_node, t_node->children[0], t_node->p_k, t_node->p_j, useless_point_seq);
                        // adjacent i,j,p
                        t_node->children[0]->adjacent[1] = t_node->children[1];
                        // adjacent k,l,p
                        if(t_adjacent != NULL)
                            t_node->children[0]->adjacent[2] = t_adjacent->children[0];

                        // triangle i,j,p:
                        // adjacent by parent
                        t_node->children[1]->adjacent[0] = adjacent_triangle_by_parent(t_node, t_node->children[1], t_node->p_i, t_node->p_j, useless_point_seq);
                        // adjacent k,j,p
                        t_node->children[1]->adjacent[1] = t_node->children[0];
                        // adjacent i,l,p
                        if(t_adjacent != NULL)
                            t_node->children[1]->adjacent[2] = t_adjacent->children[1];
                        
                        if(t_adjacent != NULL)
                        {
                            // triangle k,l,p:
                            // adjacent by parent
                            t_adjacent->children[0]->adjacent[0] = adjacent_triangle_by_parent(t_adjacent, t_adjacent->children[0], point_seq[1], point_seq[0], useless_point_seq);
                            // adjacent i,l,p
                            t_adjacent->children[0]->adjacent[1] = t_adjacent->children[1];
                            // adjacent k,j,p
                            t_adjacent->children[0]->adjacent[2] = t_node->children[0];
                            //
                            // triangle i,l,p:
                            // adjacent by parent
                            t_adjacent->children[1]->adjacent[0] = adjacent_triangle_by_parent(t_adjacent, t_adjacent->children[1], point_seq[2], point_seq[0], useless_point_seq);
                            // adjacent k,l,p
                            t_adjacent->children[1]->adjacent[1] = t_adjacent->children[0];
                            // adjacent i,j,p
                            t_adjacent->children[1]->adjacent[2] = t_node->children[1];
                        }

                        // legalize the edges of the divided triangle
                        legalize_edge(p[p_index], t_node->p_k, point_seq[0], t_adjacent->children[0]->adjacent[0], t_adjacent->children[0]);
                        legalize_edge(p[p_index], point_seq[0], t_node->p_i, t_adjacent->children[1]->adjacent[0], t_adjacent->children[1]);
                        legalize_edge(p[p_index], t_node->p_i, t_node->p_j, t_node->children[1]->adjacent[0], t_node->children[1]);
                        legalize_edge(p[p_index], t_node->p_j, t_node->p_k, t_node->children[0]->adjacent[0], t_node->children[0]);
                    }

                    // no need to verify other triangles in the queue, empty the queue
                    tq = queue<struct TRIANGLE *>();
                }
            }
        }
    }

    // perform a DFS on the tree to search leaf nodes which are Delaunay triangles
    define_delaunay_triangles_from_tree(triangles_tree, triangles);

    return(triangles);
}

// save x,y coordinates of the points of each triangle in file
void save_delaunay_triangles(vector <struct TRIANGLE *> &t)
{
    int i;
    FILE *triangles_file;
    double x, y;
    
    // open file
    triangles_file = fopen("triangles.txt", "w");
    
    // process triangles
    for(i = 0;i < t.size();i++)
    {
        // save points in file
        x = t[i]->p_i.x;
        y = t[i]->p_i.y;
        fprintf(triangles_file, "%.6f %.6f\n", x, y);
        x = t[i]->p_j.x;
        y = t[i]->p_j.y;
        fprintf(triangles_file, "%.6f %.6f\n", x, y);
        x = t[i]->p_k.x;
        y = t[i]->p_k.y;
        fprintf(triangles_file, "%.6f %.6f\n", x, y);
        fprintf(triangles_file, "END\n");
    }
    
    // close file
    fclose(triangles_file);

    return;
}
