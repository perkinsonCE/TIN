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

#include <iostream>

#include "globals.h"
#include "functions.h"

// Delaunay triangles vector
vector <struct TRIANGLE *> delaunay_triangles;

int main(int argc, char **argv)
{
    // define vector of points
    cout << "Read input points..." << endl;
    read_points();

    // define initial triangle that encloses all points and
    // stores this triangle as a node of the triangles tree for the search
    // of the triangle that contains each point inserted into the triangulation
    cout << "Define enclosing triangle of input points..." << endl;
    define_enclosing_triangle(points);

    // compute the Delaunay triangulation
    cout << "Compute the Delaunay triangulation..." << endl;
    delaunay_triangles = delaunay_triangulation(points);

    cout << "Save Delaunay triangles edges in file..." << endl;
    // save Delaunay triangles edges
    save_delaunay_triangles(delaunay_triangles);

    return(0);
}
