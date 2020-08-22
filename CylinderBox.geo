boxLength = 2.5;
boxWidth = 0.41;
boxHeight = 0.41;

mSize=boxWidth/4.0;
h_refine = mSize/2.0; //refinement size



//Define Cylinder

//SetFactory("OpenCASCADE");
//Cylinder(93) = {0.1, 0.3, 4.9, 1, 0, 0, 0.5, 2*Pi};

tdisk = 0.0;
ddisk = 0.05; //radius
xoff = 0.5;
yoff = 0.2;

Point(101) = {xoff-ddisk,yoff,0.0,mSize};
Point(102) = {xoff,yoff-ddisk,0.0,mSize};
Point(103) = {xoff+ddisk,yoff,0.0,mSize};
Point(104) = {xoff,yoff+ddisk,0.0,mSize};
Point(110) = {xoff,yoff,0.0,mSize};

Point(105) = {xoff-ddisk,yoff,boxHeight,mSize};
Point(106) = {xoff,yoff-ddisk,boxHeight,mSize};
Point(107) = {xoff+ddisk,yoff,boxHeight,mSize};
Point(108) = {xoff,yoff+ddisk,boxHeight,mSize};
Point(111) = {xoff,yoff,boxHeight,mSize};

Circle(201) = {101,110,102};
Circle(202) = {102,110,103};
Circle(203) = {103,110,104};
Circle(204) = {104,110,101};
Line(105) = {101,105};
Line(106) = {102,106};
Line(107) = {103,107};
Line(108) = {104,108};
Circle(209) = {105,111,106};
Circle(210) = {106,111,107};
Circle(211) = {107,111,108};
Circle(212) = {108,111,105};
 
Line Loop(300) = {201,202,203,204};    //Plane Surface(123) = {122};
Line Loop(124) = {201,106,-209,-105};  Ruled  Surface(140) = {124};
Line Loop(126) = {202,107,-210,-106};  Ruled  Surface(141) = {126};
Line Loop(128) = {203,108,-211,-107};  Ruled  Surface(142) = {128};
Line Loop(130) = {204,105,-212,-108};  Ruled  Surface(143) = {130};
Line Loop(301) = {209,210,211,212};    //Plane Surface(133) = {132};
 
//Surface Loop(140) = {125,127,129,131} ;
Surface Loop(140) = {140,141,142,143} ;

//mSize=boxLength/4.0;

Point(58) = {0.0, 0.0, 0, mSize};
Point(5) = {0.0, 0.0, boxHeight,mSize};
Point(10) = {boxLength, 0.0, boxHeight,mSize};
Point(56) = {boxLength, 0.0, 0,mSize};

Point(60) = {0.0, boxWidth, 0,mSize};
Point(2) = {0.0, boxWidth, boxHeight,mSize};
Point(15) = {boxLength, boxWidth, boxHeight,mSize};
Point(54) = {boxLength, boxWidth, 0,mSize};

Line(50) = {58, 56};
Line(48) = {56, 54};
Line(46) = {54, 60};
Line(52) = {60, 58};

Line(73) = {58, 5};
Line(72) = {56, 10};
Line(71) = {54, 15};
Line(74) = {60, 2};

Line(11) = {5, 10};
Line(16) = {10, 15};
Line(20) = {15, 2};
Line(6) = {2, 5};

Line Loop(42) = {-52, -46, -48, -50}; //bottom z-
Line Loop(76) = {46,74,-20,-71}; //back y+
Line Loop(24) = {16,20,6,11}; //top z+
Line Loop(80) = {50,72,-11,-73}; //front y-
Line Loop(82) = {-6,-74,52,73}; // left x-
Line Loop(78) = {48, 71, -16, -72}; // right x+

Plane Surface(1) = {42,300}; //bottom 1
Plane Surface(5) = {80}; //front 
Plane Surface(4) = {78}; //right
Plane Surface(3) = {76}; //back 
Plane Surface(2) = {82}; //left
Plane Surface(6) = {24,301}; //top

Surface Loop(51) = {1,5,4,3,2,6};


Volume(92) = {51,140};
//+
Coherence;

//Specifying the physical entities removes all mesh elements that are not part of the physical entities, so it removes unused points like that used to define the sphere

Physical Point(150) = {58,5,10,56,60,2,15,54};
//Physical Point(151) = {101,102,103,104,105,106,107,108};
//Physical Point(150) = {58,5,10,56,60,2,15,54,101,102,103,104,105,106,107,108};
Physical Line(160)  = {50,48,46,52,73,72,71,74,11,16,20,6,201,202,203,204,209,210,211,212};
//Physical Line(161) = {105,106,107,108}
//Physical Surface (170) = {1,5,4,3,2,6,125,127,129,131};
Physical Surface (170) = {1,5,4,3,2,6,140,141,142,143};
Physical Volume(1000) = {92};

Coherence;
//+

