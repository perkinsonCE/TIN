# TIN

This C++ implementation of the Delaunay triangulation is available with a simple makefile.
After extracting the source code, you just need to enter the command "make TIN" (Unix systems).

The program reads the "points.txt" file that contains all points to be triangulated in the following pattern:

n<br>
1 X1 Y1 Z1<br>
2 X2 Y2 Z2<br>
...<br>
n Xn Yn Zn

where the first line contains the number of points "n" and the next "n" lines contain all points specified by
an index and three plane coordinates (X,Y,Z).

The execution of the TIN program creates the "triangles.txt" file with the plane coordinates (X,Y) of all Delaunay
triangles. All triangles are separated by the word "END".
