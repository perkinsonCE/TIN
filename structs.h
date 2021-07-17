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

struct POINT
{
    // point coordinates
    double x, y, z;
    // point memory address
    struct POINT *address;
    // point index in vector
    int point_index;
    // flags
    bool processed, visited;
};

struct TRIANGLE
{
    // points of the triangle
    struct POINT p_i, p_j, p_k;
    // number of children and parents of the triangle
    int number_of_children, number_of_parents;
    // references to the children, adjacent and parents of
    // the triangle
    struct TRIANGLE **children, *adjacent[3], *parents[2];
    // triangle index in vector
    int triangle_index;
    // circumcenter coordinates of the circunference defined
    // from the 3 points of the triangle
    double xc, yc;
    // circunference radius
    double r;
    // flags
    bool processed, visited;
};
