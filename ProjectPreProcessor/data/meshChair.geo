a = 5;
lc = 0.1;

Point(1) = {0, 0, 0, lc};

Point(3) = { a*0.2, 0, 0, lc};
Point(4) = {a*0.2, a*-0.8, 0, lc};
Point(5) = {a*1,a*-1,0,lc}; 
Point(6) = {a*1, a*-0.8, 0, lc};
Point(7) = {a*0.1, 0, 0, lc};
Point (8) = {0,a*-0.8,0,lc};
Point (9) = {a*0.2,a*-1,0,lc};
Point (10) = {2, -5, 0, lc};
Point (11) = {3.25, -5, 0, lc};
Point(12) = {2,-7,0,lc};
Point(13) = {3.25,-7,0,lc};

Point(19) = {1, -3.6, 0, lc};
Point(20) = {1.4, -4.0269, 0, lc};




Line(1) = {12,13};
Line(2) = {13,11};
Line(3) = {11,5};
Line(4) = {5,6};
Line(5) = {6,20};
Bezier(6) = {19, 4, 4, 20};
Line(7) = {19,3};
Circle (8) = {3,7,1};
Line(9) = {1,8};
Circle (10) = {8,4,9};
Line(11) = {9,10};
Line(12) = {10,12};









Curve Loop(1) = {1:5,-6,7:12};

Plane Surface(1000) = {1};
Recombine Surface (1000);