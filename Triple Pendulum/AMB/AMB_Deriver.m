% Variables needed for AMB
syms th1 th2 th3 dth1 dth2 dth3 ddth1 ddth2 ddth3 ;

% constants needed
syms l1 l2 l3 m1 m2 m3 I1 I2 I3 g d1 d2 d3 ;

% unit vectors:
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
   
% Force of gravity on each link:

Fg1 = -m1 * g * j ;
Fg2 = -m2 * g * j ;
Fg3 = -m3 * g * j ;

% Angular velocity vectors:

th1Dot = dth1 * k ;
th2Dot = dth2 * k ;
th3Dot = dth3 * k ;

% Angular acceleration vectors:

th1Ddot = ddth1 * k ;
th2Ddot = ddth2 * k ;
th3Ddot = ddth3 * k ;

% Kinematics:

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

% COM kinematics:

% Link 1:
r1_com = r1_o/ l1 * d1 ;
v1_com = v1_o/ l1 * d1 ;
a1_com = a1_o/ l1 * d1 ;

% Link 2:

r2_com = r1_o + Ie2\ r2 / l2 * d2 ;
v2_com = v1_o + Ie2\ dr2 / l2 * d2 ;
a2_com = a1_o + Ie2\ ddr2 / l2 * d2 ;

% Link 3:

r3_com = r2_o + Ie3\ r3 / l3 * d3 ;
v3_com = v2_o + Ie3\ dr3 / l3 * d3 ;
a3_com = a2_o + Ie3\ ddr3 / l3 * d3 ;

% Begin angular momentuum balance:

% AMB 1:
H1 = m1 * cross( r1_com, a1_com ) + I1 * th1Ddot ;
H2 = m2 * cross( r2_com, a2_com ) + I2 * th2Ddot ;
H3 = m3 * cross( r3_com, a3_com ) + I3 * th3Ddot ;

H_o = H1 + H2 + H3 ;
M_o = cross( r1_com, Fg1 ) + cross( r2_com, Fg2 ) + cross( r3_com, Fg3 ) ;

AMB1 = simplify( H_o == M_o ) ;

% AMB 2:

H2_1 = m2 * cross( r2_com - r1_o, a2_com ) + I2 * th2Ddot ;
H3_1 = m3 * cross( r3_com - r1_o, a3_com ) + I3 * th3Ddot ;

H_1 = H2_1 + H3_1 ;
M_1 = cross( r2_com - r1_o, Fg2 ) + cross( r3_com - r1_o, Fg3 ) ;

r23g_o = ( m2 * r2_com + m3 * r3_com ) / ( m2 + m3) ;
cor1 = 0 ;

AMB2 = simplify( H_1 == M_1 + cor1 ) ;

% AMB 3:
H3_2 = m3 * cross( r3_com - r2_o, a3_com ) + I3 * th3Ddot ;

M3_2 = cross( r3_com - r2_o, Fg3 ) ;

cor2 = 0 ;

AMB3 = simplify( H3_2 == M3_2 + cor2 ) ;

% Solve for the equations of motion
EOM = solve( [ AMB1(3), AMB2(3), AMB3(3)], [ddth1, ddth2, ddth3] ) ;

% Extract equations of motion from the structure
th1DotDot = simplify( EOM.ddth1 ) ;
th2DotDot = simplify( EOM.ddth2 ) ;
th3DotDot = simplify( EOM.ddth3 ) ;

% Create matlab function for each equation
matlabFunction( th1DotDot, 'File', 'the1DotDot' ) ;
matlabFunction( th2DotDot, 'File', 'the2DotDot' ) ;
matlabFunction( th3DotDot, 'File', 'the3DotDot' ) ;