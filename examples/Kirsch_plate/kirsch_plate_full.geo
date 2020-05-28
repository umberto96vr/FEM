// Mesh size
h = 0.34;

// Constants
rd = 2.5;     // Hole radius.
lx = 20;    // Plate width.
ly = 20;    // Plate height.

// Points definition
Point(1) = {      0,       0, 0, h};
Point(2) = {     rd,       0, 0, h};
Point(3) = {      0,      rd, 0, h};
Point(4) = {    -rd,       0, 0, h};
Point(5) = {      0,     -rd, 0, h};
Point(6) = { 0.5*lx,  0.5*ly, 0, h};
Point(7) = {-0.5*lx,  0.5*ly, 0, h};
Point(8) = {-0.5*lx, -0.5*ly, 0, h};
Point(9) = { 0.5*lx, -0.5*ly, 0, h};

// Line definition
Circle(23) = {2,1,3};   //62 to have counter clock-wise direction
Circle(34) = {3,1,4};
Circle(45) = {4,1,5};
Circle(52) = {5,1,2};
Line(78)   = {7,8};
Line(89)   = {8,9};
Line(96)   = {9,6};
Line(67)   = {6,7};

// Line loop
Line Loop(1) = {78, 89, 96, 67};
Line Loop(2) = {23, 34, 45, 52};

// Surface generation
Plane Surface(1) = {1,2};

// Physical surface
Physical Surface("piastra", 100) = {1};

// Physical lines
//Physical Curve("EnforcedDisplacement",   101) = {78, 89, 96, 67, 23, 34, 45, 52};
Physical Curve("EnforcedDisplacement",   101) = {78, 89, 96, 67};