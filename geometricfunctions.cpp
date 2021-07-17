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

#include <cmath>

#include "globals.h"
#include "functions.h"

// flags to indicate if point is on an edge
bool same_slope_ij, same_slope_jk, same_slope_ki;

// calculates euclidian distance between 2 points in the plane
double distance_between(struct POINT p, struct POINT q)
{
    double dist, dx, dy;
    
    dx = p.x - q.x;
    dy = p.y - q.y;
    dist = sqrt(dx*dx+dy*dy);
    
    return(dist);
}

// calculates distance between point p and line Ax+By+C with coefficients A,B,C as |A*x+B*y+C|/SQRT(A^2+B^2)
double point_line_distance(double A, double B, double C, struct POINT p)
{
    double distance;

    distance = abs(A*(p.x) + B*(p.y) + C)/sqrt(A*A + B*B);
    
    return(distance);
}

// verifies if point p is on the line that passes through points p1,p2
//
// uses point-line distance formula by calculating the coefficients A,B,C of the line Ax+By+C=0 from points p1,p2
bool same_slope(struct POINT p, struct POINT p1, struct POINT p2)
{
    double A, B, C, distance;
    
    // coefficients of line p1,p2
    A = p1.y - p2.y;
    B = p2.x - p1.x;
    C = (p1.x - p2.x)*p1.y + (p2.y - p1.y)*p1.x;
    
    // distance between point p and line p1,p2
    distance = point_line_distance(A, B, C, p);

    // if distance is less than an error, then point is on the line
    if(distance < epsEdge)
        return(true);
    else
        return(false);
}

// calculates intersection point between 2 segments p1,p2 and p3,p4
//
// cases of vertical/horizontal lines are treated separately
//
// coordinates of intersection point are the solution of system of 2 equations
struct POINT segment_intersection(struct POINT p1, struct POINT p2, struct POINT p3, struct POINT p4)
{
    struct POINT p;
    double xi, yi, m12, m34;

    // parallel lines
    if(same_slope(p1, p3, p4) == true && same_slope(p2, p3, p4) == true)
    {
        p.x = UNDEF;
        p.y = UNDEF;
        
        return(p);
    }
    
    // segment p1,p2 is vertical
    if(abs(p1.x - p2.x) < epsVertex)
    {
        // segment p3,p4 is vertical, so both are parallel, there is no intersection point
        if(abs(p3.x - p4.x) < epsVertex)
        {
            p.x = UNDEF;
            p.y = UNDEF;
            
            return(p);
        }

        // segment p3,p4 is horizontal
        if(abs(p3.y - p4.y) < epsVertex)
        {
            p.x = p1.x;
            p.y = p3.y;
            
            if(((-epsVertex + p3.x <= p.x && p.x <= p4.x + epsVertex) || (-epsVertex + p4.x <= p.x && p.x <= p3.x + epsVertex)) &&
               ((-epsVertex + p1.y <= p.y && p.y <= p2.y + epsVertex) || (-epsVertex + p2.y <= p.y && p.y <= p1.y + epsVertex)))
                return(p);
            
            p.x = UNDEF;
            p.y = UNDEF;
            
            return(p);
        }

        // segment p3,p4 is not vertical/horizontal
        xi = p1.x;
        m34 = (p4.y - p3.y)/(p4.x - p3.x);
        yi = p3.y + m34*(xi - p3.x);
        
        p.x = xi;
        p.y = yi;
        
        if(((-epsVertex + p3.x <= p.x && p.x <= p4.x + epsVertex) || (-epsVertex + p4.x <= p.x && p.x <= p3.x + epsVertex)) &&
           ((-epsVertex + p3.y <= p.y && p.y <= p4.y + epsVertex) || (-epsVertex + p4.y <= p.y && p.y <= p3.y + epsVertex)) &&
           ((-epsVertex + p1.y <= p.y && p.y <= p2.y + epsVertex) || (-epsVertex + p2.y <= p.y && p.y <= p1.y + epsVertex)))
            return(p);
            
        p.x = UNDEF;
        p.y = UNDEF;

        return(p);
    }

    // segment p1,p2 is horizontal
    if(abs(p1.y - p2.y) < epsVertex)
    {
        // segment p3,p4 is horizontal, so both are parallel, there is no intersection point
        if(abs(p3.y - p4.y) < epsVertex)
        {
            p.x = UNDEF;
            p.y = UNDEF;

            return(p);
        }

        // segment p3,p4 is vertical
        if(abs(p3.x - p4.x) < epsVertex)
        {
            p.x = p3.x;
            p.y = p1.y;
            
            if(((-epsVertex + p1.x <= p.x && p.x <= p2.x + epsVertex) || (-epsVertex + p2.x <= p.x && p.x <= p1.x + epsVertex)) &&
               ((-epsVertex + p3.y <= p.y && p.y <= p4.y + epsVertex) || (-epsVertex + p4.y <= p.y && p.y <= p3.y + epsVertex)))
                return(p);
            
            p.x = UNDEF;
            p.y = UNDEF;
            
            return(p);
        }

        // segment p3,p4 is not vertical/horizontal
        yi = p1.y;
        m34 = (p4.y - p3.y)/(p4.x - p3.x);
        xi = (yi - p3.y)/m34 + p3.x;
        
        p.x = xi;
        p.y = yi;

        if(((-epsVertex + p3.x <= p.x && p.x <= p4.x + epsVertex) || (-epsVertex + p4.x <= p.x && p.x <= p3.x + epsVertex)) &&
           ((-epsVertex + p3.y <= p.y && p.y <= p4.y + epsVertex) || (-epsVertex + p4.y <= p.y && p.y <= p3.y + epsVertex)) &&
           ((-epsVertex + p1.x <= p.x && p.x <= p2.x + epsVertex) || (-epsVertex + p2.x <= p.x && p.x <= p1.x + epsVertex)))
            return(p);
            
        p.x = UNDEF;
        p.y = UNDEF;

        return(p);
    }

    // segment p1,p2 is not vertical/horizontal
    //
    // segment p3,p4 is vertical
    if(abs(p3.x - p4.x) < epsVertex)
    {
        xi = p3.x;
        m12 = (p2.y - p1.y)/(p2.x - p1.x);
        yi = p1.y + m12*(xi - p1.x);
    
        p.x = xi;
        p.y = yi;
    
        if(((-epsVertex + p1.x <= p.x && p.x <= p2.x + epsVertex) || (-epsVertex + p2.x <= p.x && p.x <= p1.x + epsVertex)) &&
           ((-epsVertex + p1.y <= p.y && p.y <= p2.y + epsVertex) || (-epsVertex + p2.y <= p.y && p.y <= p1.y + epsVertex)) &&
           ((-epsVertex + p3.y <= p.y && p.y <= p4.y + epsVertex) || (-epsVertex + p4.y <= p.y && p.y <= p3.y + epsVertex)))
            return(p);
        
        p.x = UNDEF;
        p.y = UNDEF;

        return(p);
    }
    // segment p3,p4 is horizontal
    if(abs(p3.y - p4.y) < epsVertex)
    {
        yi = p3.y;
        m12 = (p2.y - p1.y)/(p2.x - p1.x);
        xi = (yi - p1.y)/m12 + p1.x;
        
        p.x = xi;
        p.y = yi;
		
        if(((-epsVertex + p1.x <= p.x && p.x <= p2.x + epsVertex) || (-epsVertex + p2.x <= p.x && p.x <= p1.x + epsVertex)) &&
           ((-epsVertex + p1.y <= p.y && p.y <= p2.y + epsVertex) || (-epsVertex + p2.y <= p.y && p.y <= p1.y + epsVertex)) &&
           ((-epsVertex + p3.x <= p.x && p.x <= p4.x + epsVertex) || (-epsVertex + p4.x <= p.x && p.x <= p3.x + epsVertex)))
            return(p);
        
        p.x = UNDEF;
        p.y = UNDEF;
        
        return(p);
    }

    // both segments are not horizontal/vertical
    m12 = (p2.y - p1.y)/(p2.x - p1.x);
    m34 = (p4.y - p3.y)/(p4.x - p3.x);
    xi = (p3.y - p1.y - m34*p3.x + m12*p1.x)/(m12 - m34);
    yi = p1.y + m12*(xi - p1.x);
    
    p.x = xi;
    p.y = yi;
    
    if(((-epsVertex + p1.x <= p.x && p.x <= p2.x + epsVertex) || (-epsVertex + p2.x <= p.x && p.x <= p1.x + epsVertex)) &&
       ((-epsVertex + p1.y <= p.y && p.y <= p2.y + epsVertex) || (-epsVertex + p2.y <= p.y && p.y <= p1.y + epsVertex)) &&
       ((-epsVertex + p3.x <= p.x && p.x <= p4.x + epsVertex) || (-epsVertex + p4.x <= p.x && p.x <= p3.x + epsVertex)) &&
       ((-epsVertex + p3.y <= p.y && p.y <= p4.y + epsVertex) || (-epsVertex + p4.y <= p.y && p.y <= p3.y + epsVertex)))
        return(p);
        
    p.x = UNDEF;
    p.y = UNDEF;

    return(p);
}

// verifies if point is inside triangle, on an edge or outside triangle
//
// considers a horizontal line from the point going to the right and checks how many and what kind of intersections with edges occur
//
// return values
// 0 - point is outside triangle
// 1 - point is on an edge
// 2 - point is inside triangle
int point_in_triangle(struct POINT p, struct TRIANGLE t)
{
    struct POINT pmax, intersection_point;
    double max_x;
    int edge_intersections;
    bool up, down;
    
    // if point is on the edge i,j
    if(same_slope(p, t.p_i, t.p_j) == true)
    {
        same_slope_ij = true;
        same_slope_jk = same_slope_ki = false;

        // if coordinates of the point are in the range of the coordinates of the edge points
        if(((-epsVertex + t.p_i.x <= p.x && p.x <= t.p_j.x + epsVertex) || (-epsVertex + t.p_j.x <= p.x && p.x <= t.p_i.x + epsVertex)) &&
           ((-epsVertex + t.p_i.y <= p.y && p.y <= t.p_j.y + epsVertex) || (-epsVertex + t.p_j.y <= p.y && p.y <= t.p_i.y + epsVertex)))
            return(1);
        else
            return(0);
    }
    
    // if point is on the edge j,k
    if(same_slope(p, t.p_j, t.p_k) == true)
    {
        same_slope_jk = true;
        same_slope_ki = same_slope_ij = false;

        // if coordinates of the point are in the range of the coordinates of the edge points
        if(((-epsVertex + t.p_j.x <= p.x && p.x <= t.p_k.x + epsVertex) || (-epsVertex + t.p_k.x <= p.x && p.x <= t.p_j.x + epsVertex)) &&
           ((-epsVertex + t.p_j.y <= p.y && p.y <= t.p_k.y + epsVertex) || (-epsVertex + t.p_k.y <= p.y && p.y <= t.p_j.y + epsVertex)))
            return(1);
        else
            return(0);
    }

    // if point is on the edge k,i
    if(same_slope(p, t.p_k, t.p_i) == true)
    {
        same_slope_ki = true;
        same_slope_ij = same_slope_jk = false;

        // if coordinates of the point are in the range of the coordinates of the edge points
        if(((-epsVertex + t.p_k.x <= p.x && p.x <= t.p_i.x + epsVertex) || (-epsVertex + t.p_i.x <= p.x && p.x <= t.p_k.x + epsVertex)) &&
           ((-epsVertex + t.p_k.y <= p.y && p.y <= t.p_i.y + epsVertex) || (-epsVertex + t.p_i.y <= p.y && p.y <= t.p_k.y + epsVertex)))
            return(1);
        else
            return(0);
    }
    
    // defines maximum x coordinate from points i,j,k of the triangle
    max_x = -1.;
    if(t.p_i.x > max_x)
        max_x = t.p_i.x;
    if(t.p_j.x > max_x)
        max_x = t.p_j.x;
    if(t.p_k.x > max_x)
        max_x = t.p_k.x;

    // coordinates of the outside point pmax
    pmax.x = max_x + INC;
    pmax.y = p.y;
                
    // verifies if point is inside triangle by the intersection between horizontal line from point p with the triangle edges
    //
    // initialize number of intersections and direction of edge points
    edge_intersections = 0;
    up = down = false;
    // intersection between segments p,pmax and i,j
    intersection_point = segment_intersection(p, pmax, t.p_i, t.p_j);
    // if intersection exists
    if(point_not_undef(intersection_point) == true)
    {
        // increment number of intersections
        edge_intersections++;
         
        // if intersection is a vertex
        if(points_are_equal(t.p_i, intersection_point) == true)
        {
            // if intersected vertex is up or down on the edge
            if(t.p_i.y < t.p_j.y)
                down = true;
            else
                up = true;
        }
        else
        // if intersection is a vertex
        if(points_are_equal(t.p_j, intersection_point) == true)
        {
            // if intersected vertex is up or down on the edge
            if(t.p_j.y < t.p_i.y)
                down = true;
            else
                up = true;
        }
    }
        
    // intersection between segments p,pmax and j,k
    intersection_point = segment_intersection(p, pmax, t.p_j, t.p_k);
    // if intersection exists
    if(point_not_undef(intersection_point) == true)
    {
        // increment number of intersections
        edge_intersections++;
         
        // if intersection is a vertex
        if(points_are_equal(t.p_j, intersection_point) == true)
        {
            // if intersected vertex is up or down on the edge
            if(t.p_j.y < t.p_k.y)
                down = true;
            else
                up = true;
        }
        else
        // if intersection is a vertex
        if(points_are_equal(t.p_k, intersection_point) == true)
        {
            // if intersected vertex is up or down on the edge
            if(t.p_k.y < t.p_j.y)
                down = true;
            else
                up = true;
        }
    }

    // intersection between segments p,pmax and k,i
    intersection_point = segment_intersection(p, pmax, t.p_k, t.p_i);
    // if intersection exists
    if(point_not_undef(intersection_point) == true)
    {
        // increment number of intersections
        edge_intersections++;
         
        // if intersection is a vertex
        if(points_are_equal(t.p_k, intersection_point) == true)
        {
            // if intersected vertex is up or down on the edge
            if(t.p_k.y < t.p_i.y)
                down = true;
            else
                up = true;
        }
        else
        // if intersection is a vertex
        if(points_are_equal(t.p_i, intersection_point) == true)
        {
            // if intersected vertex is up or down on the edge
            if(t.p_i.y < t.p_k.y)
                down = true;
            else
                up = true;
        }
    }
    
    // if number of intersections with edges is either 1 or 2 and it happened with a vertex forming a right "beak", then point is inside triangle
    if((edge_intersections == 1) || (edge_intersections == 2 && down == true && up == true))
        return(2);
    
    return(0);
}

// calculates the center and radius of a circunference passing through points p,p1,p2
void circle_by_3_points(struct POINT p, struct POINT p1, struct POINT p2, double *xc, double *yc, double *r)
{
    double dx, dy;
    double m1, m2, xm1, ym1, xm2, ym2;

    // verifies if line p1,p2 is vertical and p2,p is horizontal, then circumcenter is defined by mean point
    if(abs(p1.x - p2.x) < epsVertex && abs(p2.y - p.y) < epsVertex)
    {
        // calculates circumcenter coordinates
        *xc = (p2.x + p.x)/2;
        *yc = (p1.y + p2.y)/2;
    }
    else
    // verifies if line p2,p is vertical and p,p1 is horizontal, then circumcenter is defined by mean point
    if(abs(p2.x - p.x) < epsVertex && abs(p.y - p1.y) < epsVertex)
    {
        // calculates circumcenter coordinates
        *xc = (p.x + p1.x)/2;
        *yc = (p2.y + p.y)/2;
    }
    else
    // verifies if line p,p1 is vertical and p1,p2 is horizontal, then circumcenter is defined by mean point
    if(abs(p.x - p1.x) < epsVertex && abs(p1.y - p2.y) < epsVertex)
    {
        // calculates circumcenter coordinates
        *xc = (p1.x + p2.x)/2;
        *yc = (p.y + p1.y)/2;
    }
    else
    // verifies if line p1,p2 is horizontal and p2,p is vertical, then circumcenter is defined by mean point
    if(abs(p1.y - p2.y) < epsVertex && abs(p2.x - p.x) < epsVertex)
    {
        // calculates circumcenter coordinates
        *xc = (p1.x + p2.x)/2;
        *yc = (p2.y + p.y)/2;
    }
    else
    // verifies if line p2,p is horizontal and p,p1 is vertical, then circumcenter is defined by mean point
    if(abs(p2.y - p.y) < epsVertex && abs(p.x - p1.x) < epsVertex)
    {
        // calculates circumcenter coordinates
        *xc = (p2.x + p.x)/2;
        *yc = (p.y + p1.y)/2;
    }
    else
    // verifies if line p,p1 is horizontal and p1,p2 is vertical, then circumcenter is defined by mean point
    if(abs(p.y - p1.y) < epsVertex && abs(p1.x - p2.x) < epsVertex)
    {
        // calculates circumcenter coordinates
        *xc = (p.x + p1.x)/2;
        *yc = (p1.y + p2.y)/2;
    }
    else
    // verifies if line p1,p2 is vertical, then circumcenter is defined by constant line
    if(abs(p1.x - p2.x) < epsVertex)
    {
        // calculates circumcenter coordinates
        // vertical line, calculates mean y and x by intersection
        *yc = (p1.y + p2.y)/2;
        xm2 = (p2.x + p.x)/2;
        ym2 = (p2.y + p.y)/2;
        m2 = (p.y - p2.y)/(p.x - p2.x);
        *xc = -(*yc - ym2)*m2 + xm2;
    }
    else
    // verifies if line p2,p is vertical, then circumcenter is defined by constant line
    if(abs(p2.x - p.x) < epsVertex)
    {
        // calculates circumcenter coordinates
        // vertical line, calculates mean y and x by intersection
        *yc = (p2.y + p.y)/2;
        xm1 = (p1.x + p2.x)/2;
        ym1 = (p1.y + p2.y)/2;
        m1 = (p2.y - p1.y)/(p2.x - p1.x);
        *xc = -(*yc - ym1)*m1 + xm1;
    }
    else
    // verifies if line p,p1 is vertical, then circumcenter is defined by constant line
    if(abs(p.x - p1.x) < epsVertex)
    {
        // calculates circumcenter coordinates
        // vertical line, calculates mean y and x by intersection
        *yc = (p.y + p1.y)/2;
        xm1 = (p1.x + p2.x)/2;
        ym1 = (p1.y + p2.y)/2;
        m1 = (p2.y - p1.y)/(p2.x - p1.x);
        *xc = -(*yc - ym1)*m1 + xm1;
    }
    else
    // verifies if line p1,p2 is horizontal, then circumcenter is defined by constant line
    if(abs(p1.y - p2.y) < epsVertex)
    {
        // calculates circumcenter coordinates
        // vertical line, calculates mean x and y by intersection
        *xc = (p1.x + p2.x)/2;
        xm2 = (p2.x + p.x)/2;
        ym2 = (p2.y + p.y)/2;
        m2 = (p.y - p2.y)/(p.x - p2.x);
        *yc = -(*xc - xm2)/m2 + ym2;
    }
    else
    // verifies if line p2,p is horizontal, then circumcenter is defined by constant line
    if(abs(p2.y - p.y) < epsVertex)
    {
        // calculates circumcenter coordinates
        // vertical line, calculates mean x and y by intersection
        *xc = (p2.x + p.x)/2;
        xm1 = (p1.x + p2.x)/2;
        ym1 = (p1.y + p2.y)/2;
        m1 = (p2.y - p1.y)/(p2.x - p1.x);
        *yc = -(*xc - xm1)/m1 + ym1;
    }
    else
    // verifies if line p,p1 is horizontal, then circumcenter is defined by constant line
    if(abs(p.y - p1.y) < epsVertex)
    {
        // calculates circumcenter coordinates
        // vertical line, calculates mean x and y by intersection
        *xc = (p.x + p1.x)/2;
        xm1 = (p1.x + p2.x)/2;
        ym1 = (p1.y + p2.y)/2;
        m1 = (p2.y - p1.y)/(p2.x - p1.x);
        *yc = -(*xc - xm1)/m1 + ym1;
    }
    else
    {
        // if both lines are not horizontal/vertical
        xm1 = (p1.x + p2.x)/2;
        ym1 = (p1.y + p2.y)/2;
        xm2 = (p2.x + p.x)/2;
        ym2 = (p2.y + p.y)/2;
        m1 = (p2.y - p1.y)/(p2.x - p1.x);
        m2 = (p.y - p2.y)/(p.x - p2.x);
        *xc = (xm2/m2 + ym2 - xm1/m1 - ym1)/(1/m2 - 1/m1);
        *yc = -(*xc - xm1)/m1 + ym1;
    }
    
    // calculates radius
    dx = p.x - *xc;
    dy = p.y - *yc;
    *r = sqrt(dx*dx + dy*dy);

    return;
}

