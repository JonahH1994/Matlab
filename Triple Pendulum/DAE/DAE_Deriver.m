% Variables needed for DAE:
syms x1 x2 x3 y1 y2 y3 dx1 dx2 dx3 dy1 dy2 dy3 ddx1 ddx2 ddx3 real ;
syms ddy1 ddy2 ddy3 th1 th2 th3 dth1 dth2 dth3 ddth1 ddth2 ddth3 real ;
syms Fx Fy F1x F2x F1y F2y ;

ddx     = [ ddx1; ddx2; ddx3 ] ;
ddy     = [ ddy1; ddy2; ddy3 ] ;
ddth    = [ ddth1; ddth2; ddth3 ] ;
allF    = [ Fx; Fy; F1x; F1y; F2x; F2y ] ;

% Constants needed:
syms l1 l2 l3 d1 d2 d3 m1 m2 m3 I1 I2 I3 g real ;

% unit vectors:

i = [ 1; 0; 0 ] ;
j = [ 0; 1; 0 ] ;
k = [ 0; 0; 1 ] ;

% Force vectors:

F   = Fx*i + Fy*j ;
F1  = F1x*i + F1y*j ;
F2  = F2x*i + F2y*j ;

% Force of gravity: 

Fg1 = -m1*g*j ;
Fg2 = -m2*g*j ; 
Fg3 = -m3*g*j ;

% Angular velocity vector:

th1Dot = dth1 * k ;
th2Dot = dth2 * k ;
th3Dot = dth3 * k ;

% Angular acceleration vectors:

th1Ddot = ddth1 * k ;
th2Ddot = ddth2 * k ;
th3Ddot = ddth3 * k ;

% Kinematics:

% Link 1:
r1_i = [ x1; y1; 0 ] ;
v1_i = [ dx1; dy1; 0 ] ;
a1_i = [ ddx1; ddy1; 0 ] ;

% Link 2:
r2_i = [ x2; y2; 0 ] ;
v2_i = [ dx2; dy2; 0 ] ;
a2_i = [ ddx2; ddy2; 0 ] ;

% Link 3:
r3_i = [ x3; y3; 0 ] ;
v3_i = [ dx3; dy3; 0 ] ;
a3_i = [ ddx3; ddy3; 0 ] ;

r1_i_com = r1_i * d1/l1 ;
a1_i_com = a1_i * d1/l1 ;

r2_i_com = r1_i + (r2_i-r1_i) * d2/l2 ;
a2_i_com = a1_i + (a2_i-a1_i) * d2/l2 ;

r3_i_com = r2_i + (r3_i-r2_i)* d3/l3 ;
a3_i_com = a2_i + (a3_i-a2_i)* d3/l3 ;

% LMB:

eqn12 = ( m1 * a1_i_com == F + F1 + Fg1 ) ;
eqn34 = ( m2 * a2_i_com == -F1 + F2 + Fg2 ) ;
eqn56 = ( m3 * a3_i_com == -F2 + Fg3 ) ;

% AMB:

H1 = I1 * th1Ddot ;
eqn7 = ( H1 == cross( (r1_i-r1_i_com), F1 ) + cross( -r1_i_com, F ) ) ;

H2 = I2 * th2Ddot ;
eqn8 = ( H2 == cross( (r2_i-r2_i_com), F2 ) + cross( -(r2_i_com-r1_i), -F1) ) ;

H3 = I3 * th3Ddot ;
eqn9 = ( H3 == cross( -(r3_i_com-r2_i), -F2 ) ) ;

% Kinematics in polar coordinates:

% Transformation matrices:

Ie1 = [ cos(th1), sin(th1), 0; ...
       -sin(th1), cos(th1), 0; ...
       0, 0, 1 ] ;

Ie2 = [ cos(th2), sin(th2), 0 ; ...
       -sin(th2), cos(th2), 0 ; ...
       0, 0, 1 ] ;
   
Ie3 = [ cos(th3), sin(th3), 0 ; ...
       -sin(th3), cos(th3), 0 ; ...
       0, 0, 1 ] ;

% Convert everything into inertial frame:

% Link 1:
r1      = [ l1; 0; 0 ] ;
dr1     = cross( th1Dot, r1 ) ;
ddr1    = jacobian( dr1, dth1 ) * ddth1 + cross( th1Dot, dr1 ) ;

% Link 2:
r2      = [ l2; 0; 0 ] ;
dr2     = cross( th2Dot, r2 ) ;
ddr2    = jacobian( dr2, dth2 ) * ddth2 + cross( th2Dot, dr2 ) ;

% Link 3:
r3      = [ l3; 0; 0 ] ;
dr3     = cross( th3Dot, r3 ) ;
ddr3    = jacobian( dr3, dth3 ) * ddth3 + cross( th3Dot, dr3 ) ;

% Link 1:

r1_o = ( Ie1\ r1 ) ;
v1_o = ( Ie1\ dr1 ) ;
a1_o = ( Ie1\ ddr1 ) ;

% Link 2:

r2_o = ( r1_o + Ie2\ r2 ) ;
v2_o = ( v1_o + Ie2\ dr2 ) ;
a2_o = ( a1_o + Ie2\ ddr2 ) ;

% Link 3:

r3_o = r2_o + Ie3\ r3 ;
v3_o = ( v2_o + Ie3\ dr3 ) ;
a3_o = ( a2_o + Ie3\ ddr3 ) ;

% Link 1:
r1_com = r1_o/ l1 * d1 ;
v1_com = v1_o/ l1 * d1 ;
a1_com = a1_o/ l1 * d1 ;

% Link 2:

r2_com = r1_o + Ie2\ r2 / l2 * d2 ;
v2_com = v1_o + Ie2\ dr2 / l2 * d2 ;
a2_com = a1_o + ( a2_o - a1_o ) * d2/ l2 ;

% Link 3:

r3_com = simplify( r2_o + Ie3\ r3 / l3 * d3 ) ;
v3_com = simplify( v2_o + Ie3\ dr3 / l3 * d3 ) ;
a3_com = a2_o + ( a3_o - a2_o ) * d3/ l3 ;

eqn10_1 = a1_com == a1_i_com ;
eqn12_3 = a2_com == a2_i_com ;
eqn14_5 = a3_com == a3_i_com ;

[ A, b] = equationsToMatrix( [eqn12(1:2); eqn34(1:2); eqn56(1:2); eqn7(3); ...
    eqn8(3); eqn9(3); eqn10_1(1:2); eqn12_3(1:2); eqn14_5(1:2)], ...
    [ ddx; ddy; ddth; allF ] ) ;

matlabFunction( A, 'File', 'AMatrix' ) ;
matlabFunction( b, 'File', 'bVector' ) ;