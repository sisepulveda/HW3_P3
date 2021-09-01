// Gmsh project created on Mon May 10 20:05:33 2021
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {20, 0, 0, 1.0};
//+
Point(3) = {20, 4, 0, 1.0};
//+
Point(4) = {0, 4, 0, 1.0};
//+
Point(5) = {2, 0, 0, 1.0};
//+
Point(6) = {18, 0, 0, 1.0};
//+
Point(7) = {18, 4, 0, 1.0};
//+
Point(8) = {2, 4, 0, 1.0};
//+
Line(1) = {1, 5};
//+
Line(2) = {5, 8};
//+
Line(3) = {8, 4};
//+
Line(4) = {4, 1};
//+
Line(5) = {5, 6};
//+
Line(6) = {6, 7};
//+
Line(7) = {7, 8};
//+
Line(8) = {6, 2};
//+
Line(9) = {2, 3};
//+
Line(10) = {3, 7};
//+
Circle(11) = {10, 2, 0, 1, 0, 2*Pi};
//+
Curve Loop(1) = {7, -2, 5, 6};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {11};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {2, 3, 4, 1};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {10, -6, 8, 9};
//+
Plane Surface(4) = {4};
//+
BooleanDifference{ Surface{1}; Delete; }{ Surface{2}; Delete; }
//+
Physical Curve("Empotrado", 1) = {4};
//+
Physical Curve("BordeNatural", 2) = {9};
//+
Physical Surface("Placa", 3) = {1};
//+
Physical Surface("Extremos", 4) = {4, 3};
//+
MeshSize {9} = 0.3;
//+
MeshSize {8, 4, 1, 5, 6, 7, 3, 2} = 1.5;
