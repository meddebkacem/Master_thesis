// Gmsh project created on Wed Apr 16 15:57:18 2025
SetFactory("OpenCASCADE");
//+
//h = 2.5; L = 10; W = 12; Ht = 4; Lt = 0.4; Wt = 0.4; W2 = 8;
h = 2; L = 20; W = 24; Ht = 8; Lt = 0.8; Wt = 0.8; W2 = 16;
//+
Point(1) = {0, 0, 0, h};
Point(2) = {L, 0, 0, h};
Point(3) = {L, Ht, 0, h/4};
Point(4) = {L+Lt,Ht,0,h/4};
Point(5) = {L+Lt,0,0, h};
Point(6) = {L+Lt+W2,0,0,h};
Point(7) = {L+Lt+W2,W,0,h};
Point(8) = {L+Lt,W,0,h};
Point(9) = {L+Lt,Ht+Wt,0,h/4};
Point(10) = {L,Ht+Wt,0,h/4};
Point(11) = {L, W, 0, h};
Point(12) = {0, W, 0, h};

ii=1;
For(1:11)
   Line(ii)={ii,ii+1};
   ii++;
EndFor
Line(12) = {12,1};
Line(13) = {3, 10};
Line(14) = {9, 4};
//+
// remove "Line(13)" & "Line(14)" to reactivate large Curve Loop
//Curve Loop(1) = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
Curve Loop(1) = {1, 2, 13, 10, 11, 12};
Curve Loop(2) = {3, -14, 9, -13};
Curve Loop(3) = {4, 5, 6, 7, 8, 14};


//+
Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};


Transfinite Line(13) = 2;
Transfinite Line(14) = 2;
Transfinite Line(3) = 2;
Transfinite Line(9) = 2;

Mesh 2;//+
Recombine Surface {1};
//+
Recombine Surface {1};
//+
Recombine Surface {1};
//+
Recombine Surface {1};



