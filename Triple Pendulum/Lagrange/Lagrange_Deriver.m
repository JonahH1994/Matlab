clear ;

% Variables needed:
syms th1 th2 th3 dth1 dth2 dth3 ddth1 ddth2 ddth3 real ;

% Constants needed:
syms l1 l2 l3 d1 d2 d3 m1 m2 m3 I1 I2 I3 g real ;

% unit vectors ;
i = [ 1; 0; 0 ] ;
j = [ 0; 1; 0 ] ;
k = [ 0; 0; 1 ] ;

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

% Kinematics in polar coordinates
  
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

% Convert everything into inertial frame:

% Link 1:
r1_o = simplify( Ie1\ r1 ) ;
v1_o = simplify( Ie1\ dr1 ) ;
a1_o = simplify( Ie1\ ddr1 ) ;

% Link 2:
r2_o = simplify( r1_o + Ie2\ r2 ) ;
v2_o = simplify( v1_o + Ie2\ dr2 ) ;
a2_o = simplify( a1_o + Ie2\ ddr2 ) ;

% Link 3:
r3_o = simplify( r2_o + Ie3\ r3 ) ;
v3_o = simplify( v2_o + Ie3\ dr3 ) ;
a3_o = simplify( a2_o + Ie3\ ddr3 ) ;

% COM equations:

% Link 1:
r1_com = r1_o/ l1 * d1 ;
v1_com = v1_o/ l1 * d1 ;
a1_com = a1_o/ l1 * d1 ;

% Link 2:

r2_com = r1_o + (r2_o-r1_o) / l2 * d2 ;
v2_com = v1_o + (v2_o-v1_o) / l2 * d2 ;
a2_com = a1_o + (a2_o-a1_o) / l2 * d2 ;

% Link 3:

r3_com = r2_o + (r3_o-r2_o) / l3 * d3 ;
v3_com = v2_o + (v3_o-v2_o) / l3 * d3 ;
a3_com = a2_o + (a3_o-a2_o) / l3 * d3 ;

% Begin deriving lagrange:

T1 = ( 0.5*m1*dot( v1_com, v1_com ) + 0.5*I1*dot( th1Dot, th1Dot ) ) ;
V1 = ( m1 * g * d1 * sin( th1 ) ) ;

T2 = ( 0.5*m2*dot( v2_com, v2_com ) + 0.5*I2*dot( th2Dot, th2Dot ) ) ;
V2 = ( m2 * g * ( l1 * sin( th1 ) + d2 * sin( th2 ) ) ) ;

T3 = ( 0.5*m3*dot( v3_com, v3_com ) + 0.5*I3*dot( th3Dot, th3Dot ) ) ;
V3 = ( m3 * g * ( l1 * sin( th1 ) + l2 * sin( th2 ) + d3 * sin( th3 ) ) ) ;

T = T1 + T2 + T3 ;
V = V1 + V2 + V3 ;

L = T - V ;

% Lagrange for th1:

diff1 = jacobian( L, dth1 ) ;
dL_dth1 = jacobian( diff1, th1 ) * dth1 + jacobian( diff1, th2 ) * dth2 + ...
    jacobian( diff1, th3 ) * dth3 + jacobian( diff1, dth1 ) * ddth1 + ...
    jacobian( diff1, dth2 ) * ddth2 + jacobian( diff1, dth3 ) * ddth3 ;

dL_dth1 = simplify( dL_dth1 ) ;

dL_th1 = jacobian( L, th1 ) ;

Q1 = 0 ;

eqn1 = simplify( dL_dth1 - dL_th1 == Q1 ) ;

% Lagrange for th2:

diff2 = jacobian( L, dth2 ) ;
dL_dth2 = jacobian( diff2, th1 ) * dth1 + jacobian( diff2, th2 ) * dth2 + ...
    jacobian( diff2, th3 ) * dth3 + jacobian( diff2, dth1 ) * ddth1 + ...
    jacobian( diff2, dth2 ) * ddth2 + jacobian( diff2, dth3 ) * ddth3 ;

dL_th2 = jacobian( L, th2 ) ;

Q2 = 0 ;

eqn2 = simplify( dL_dth2 - dL_th2 == Q2 ) ;

% Lagrange for th3:

diff3 = jacobian( L, dth3 ) ;
dL_dth3 = jacobian( diff3, th1 ) * dth1 + jacobian( diff3, th2 ) * dth2 + ...
    jacobian( diff3, th3 ) * dth3 + jacobian( diff3, dth1 ) * ddth1 + ...
    jacobian( diff3, dth2 ) * ddth2 + jacobian( diff3, dth3 ) * ddth3 ;

dL_th3 = jacobian( L, th3 ) ;

Q3 = 0 ;

eqn3 = simplify( dL_dth3 - dL_th3 == Q3 ) ;

EOM = solve( [eqn1, eqn2, eqn3], [ddth1, ddth2, ddth3] ) ;

th1DotDot = simplify( EOM.ddth1 ) ;
th2DotDot = simplify( EOM.ddth2 ) ;
th3DotDot = simplify( EOM.ddth3 ) ;

matlabFunction( th1DotDot, 'File', 'the1DotDot' ) ;
matlabFunction( th2DotDot, 'File', 'the2DotDot' ) ;
matlabFunction( th3DotDot, 'File', 'the3DotDot' ) ;