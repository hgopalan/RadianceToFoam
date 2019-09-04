lc=2.5;
Point(1)={-50,-50,0,lc};
Point(2)={50,-50,0,lc};
Point(3)={-50,50,0,lc};
Point(4)={50,50,0,lc};
lc=0.25;
Point(5)={-5,-5,0,lc};
Point(6)={5,-5,0,lc};
Point(7)={-5,5,0,lc};
Point(8)={5,5,0,lc};
//+
Line(1) = {1, 3};
//+
Line(2) = {3, 4};
//+
Line(3) = {4, 2};
//+
Line(4) = {2, 1};
//+
Line(5) = {5, 7};
//+
Line(6) = {7, 8};
//+
Line(7) = {8, 6};
//+
Line(8) = {6, 5};
//+
Extrude {0, 0, 5} {
  Line{5, 6, 7, 8};
}
//+
Extrude {0, 0, 50} {
  Line{1, 2, 3, 4};
}
//+
Line Loop(41) = {25, 29, 33, 37};
//+
Plane Surface(42) = {41};
//+
Line Loop(43) = {1, 2, 3, 4};
//+
Line Loop(44) = {5, 6, 7, 8};
//+
Plane Surface(45) = {43, 44};
//+
Line Loop(46) = {9, 13, 17, 21};
//+
Plane Surface(47) = {46};
//+
Surface Loop(48) = {28, 45, 32, 36, 40, 42, 24, 12, 16, 20, 47};
//+
Volume(49) = {48};
//+
Physical Surface("building") = {12, 24, 47, 20, 16};
//+
Physical Volume("fluid") = {49};
//+
Physical Surface("lower") = {45};
