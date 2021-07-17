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

#include "globals.h"
#include "functions.h"

// points vector
vector <struct POINT> points;

// number of points
int N;
// file of points
FILE *points_file;

void read_points(void)
{
    int i;
    int point_number;

    // open file of points
    points_file = fopen("points.txt", "r");

    // read number of points
    fscanf(points_file, "%d", &N);

    // resize vector to store all points
    points.resize(N);
    
    // read points
    for(i = 0;i < N;i++)
    {
        // read point: point number, x coordinate, y coordinate, z coordinate
        fscanf(points_file, "%d %lf %lf %lf", &point_number, &(points[i].x), &(points[i].y), &(points[i].z));
    }
    
    // close file of points
    fclose(points_file);

    // initialize information of points
    for(i = 0;i < N;i++)
    {
        // store memory address of point
        points[i].address = &(points[i]);

        // set point index
        points[i].point_index = i;
        
        // set point as not processed/visited
        points[i].processed = false;
        points[i].visited = false;
    }

    return;
}
