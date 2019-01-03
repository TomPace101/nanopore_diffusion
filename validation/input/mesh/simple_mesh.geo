//Generate mesh
//usage:
//gmsh thisfile.geo

//Mesh density control
mcarl=5.0*0.1;
mcarh=5.0*0.1;

//Geometric inputs
Lx=1.0;
Ly=1.0;
alpha1=0.5;

//Point ordinates
X1=0;
X2=Lx;
Y1=0;
Y2=(1-alpha1)*Ly;
Y3=Ly;

//Parametric locations output
Printf("# Parametric locations output") > meshmetafile;
Printf("X1: %f", X1) >> meshmetafile;
Printf("X2: %f", X2) >> meshmetafile;
Printf("Y1: %f", Y1) >> meshmetafile;
Printf("Y2: %f", Y2) >> meshmetafile;
Printf("Y3: %f", Y3) >> meshmetafile;
Printf("Xmid: %f", Lx/2) >> meshmetafile;

//Points
p11=newp; Point(p11)={X1,Y1,0,mcarl};
p12=newp; Point(p12)={X1,Y2,0,mcarh};
p13=newp; Point(p13)={X1,Y3,0,mcarl};
p21=newp; Point(p21)={X2,Y1,0,mcarl};
p22=newp; Point(p22)={X2,Y2,0,mcarh};
p23=newp; Point(p23)={X2,Y3,0,mcarl};

//Lines
L12_13=newl; Line(L12_13)={p12, p13};
L13_23=newl; Line(L13_23)={p13, p23};
L23_22=newl; Line(L23_22)={p23, p22};
L22_12=newl; Line(L22_12)={p22, p12};
L11_12=newl; Line(L11_12)={p11, p12};
L22_21=newl; Line(L22_21)={p22, p21};
L21_11=newl; Line(L21_11)={p21, p11};

//Circles

//Line Loops
Line Loop(1)={L12_13, L13_23, L23_22, L22_12};
Line Loop(2)={L11_12, -L22_12, L22_21, L21_11};

//Surfaces (each is a physical surface)
Plane Surface(1)={1};
Physical Surface(1)={1};
Plane Surface(2)={2};
Physical Surface(2)={2};

//Physical Lines
Physical Line(1213)={L12_13};
Physical Line(1323)={L13_23};
Physical Line(2223)={L23_22};
Physical Line(1222)={L22_12};
Physical Line(1112)={L11_12};
Physical Line(2122)={L22_21};
Physical Line(1121)={L21_11};

//Mesh
Mesh 2;

//Blank lines
//because gmsh gets confused without them




