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

using namespace std;

// utilityfunctions.cpp
bool point_not_undef(struct POINT p);
bool points_are_equal(struct POINT p1, struct POINT p2);
bool edges_are_equal(struct POINT p1, struct POINT p2, struct POINT p3, struct POINT p4);
struct TRIANGLE *adjacent_triangle_by_parent(struct TRIANGLE *parent, struct TRIANGLE *child, struct POINT p1, struct POINT p2, struct POINT *point_seq);
struct POINT opposite_point_by_edge(struct POINT p1, struct POINT p2, struct TRIANGLE *t);

// geometricfunctions.cpp
double distance_between(struct POINT p, struct POINT q);
double point_line_distance(double A, double B, double C, struct POINT p);
bool same_slope(struct POINT p, struct POINT p1, struct POINT p2);
struct POINT segment_intersection(struct POINT p1, struct POINT p2, struct POINT p3, struct POINT p4);
int point_in_triangle(struct POINT p, struct TRIANGLE t);
void circle_by_3_points(struct POINT p, struct POINT p1, struct POINT p2, double *xc, double *yc, double *r);

// enclosingtriangle.cpp
void define_enclosing_triangle(vector <struct POINT> &p);

// delaunaytriangulation.cpp
bool valid_delaunay_triangle(struct TRIANGLE *t);
void define_delaunay_triangles_from_tree(struct TRIANGLE *t, vector <struct TRIANGLE *> &triangles);
void legalize_edge(struct POINT p, struct POINT p1, struct POINT p2, struct TRIANGLE *t_o, struct TRIANGLE *t_p);
vector <struct TRIANGLE *> delaunay_triangulation(vector <struct POINT> &p);
void save_delaunay_triangles(vector <struct TRIANGLE *> &t);

// points.cpp
void read_points(void);
