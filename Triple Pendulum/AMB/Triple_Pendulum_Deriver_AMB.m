% free variables:
syms thet1 thet2 thet3 dthet1 dthet2 dthet3 ddthet1 ddthet2 ddthet3 ;

% constants:

syms I1 I2 I3 l1 l2 l3 m1 m2 m3 g ;

% Coordinate frame transformations:

Ie1 = [ cos(thet1), sin(thet1), 0 ; ...
    -sin(thet1), cos(thet1), 0 ; ...
    0, 0, 1 ] ;

Ie2 = [ cos(thet2), sin(thet2),0 ; ...
    -sin(thet2), cos(thet2), 0 ; ...
    0, 0, 1 ] ;

Ie3 = [ cos(thet3),  sin(thet3), 0 ; ...
    -sin(thet3), cos(thet3), 0 ; ...
    0, 0, 1 ] ;

% Unit vectors:

i = [ 1; 0; 0 ] ;
j = [ 0; 1; 0 ] ;
k = [ 0; 0; 1 ] ;

% Force of gravity: 

Fg1 = m1*g*i ;
Fg2 = m2*g*i ; 
Fg3 = m3*g*i ;

% Angular velocity vectors:

theDot1 = dthet1 * k ;
theDot2 = dthet2 * k ;
theDot3 = dthet3 * k ;

% Angular acceleration vectors:

theDdot1 = jacobian(theDot1, dthet1 ) * ddthet1 ;
theDdot2 = jacobian(theDot2, dthet2 ) * ddthet2 ;
theDdot3 = jacobian(theDot3, dthet3 ) * ddthet3 ;

% Kinematics:

r1_o = [ l1; 0; 0 ] ; % [ er1; eth1; k ] ;
r2_1 = [ l2; 0; 0 ] ; % [ er2; eth2; k ] ;
r3_2 = [ l3; 0; 0 ] ; % [ er3; eth3; k ] ;

dr1_o = cross( r1_o , theDot1) ;
dr2_1 = cross( r2_1 , theDot2) ;
dr3_2 = cross( r3_2 , theDot3) ;

ddr1_o = jacobian( dr1_o, dthet1 ) * ddthet1 + cross( dr1_o , theDot1) ;
ddr2_1 = jacobian( dr2_1, dthet2 ) * ddthet2 + cross( dr2_1 , theDot2) ;
ddr3_2 = jacobian( dr3_2, dthet3 ) * ddthet3 + cross( dr3_2 , theDot3) ;

%% Angular momentum balance:

% Everything will be in the inertial coordinate frame:

r1_1 = simplify( Ie1\ r1_o ) ;
r2_o = simplify( r1_1 + Ie2\ r2_1);
r3_o = simplify( r1_1 + Ie2\ r2_1 + Ie3\ r3_2) ;

r1_com = r1_1 / 2 ;
r2_com = (r1_1 + r2_o)/2 ;
r3_com = (r2_o + r3_o)/2 ;

dr1_com = simplify( cross( r1_1/2 , theDot1 ) ) ;
dr2_com = simplify(Ie2\ cross( r2_1/2 , theDot2 ) + dr1_com ) ;
dr3_com = simplify(Ie3\ cross( r3_2/2 , theDot3 ) + Ie2\ dr2_1 + dr1_com ) ;

% ddr1_com = simplify( Ie1\ (jacobian( dr1_o/2, dthet1 ) * ddthet1 + ...
%     cross( theDot1, dr1_o/2 ) ) ) ;

ddr1_com = simplify( Ie1\ ddr1_o/2 ) ;

% ddr2_com = simplify(Ie2\ (jacobian( dr2_1/2, dthet2 ) * ddthet2 + ...
%     cross( theDot2, dr2_1/2 )) + Ie1\ ddr1_o ) ;

ddr2_com = simplify( Ie2\ ddr2_1/2 + 2 * ddr1_com ) ;

% ddr3_com = simplify(Ie3\(jacobian( dr3_2/2, dthet3 ) * ddthet3 + ...
%     cross( theDot3, dr3_2/2 ) ) + Ie2\ddr2_1 + Ie1\ ddr1_o ) ;

ddr3_com = simplify( Ie3\ ddr3_2/2 + Ie2\ ddr2_1 + 2 * ddr1_com ) ;

% AMB of entire system about origin:

H1_O = m1*cross( r1_com, ddr1_com ) + m2*cross( r2_com, ddr2_com ) + ...
    m3*cross( r3_com, ddr3_com ) + I1*theDdot1 + I2*theDdot2 + ...
    I3*theDdot3;

Mo_1 = ( cross( r1_com, Fg1 ) + cross( r2_com, ...
    Fg2 ) + cross( r3_com, Fg3 ) ) ;

AMB1 = simplify(H1_O == Mo_1 ) ;

% AMB of first two links about joint 1:

H2_1 = m2*cross( r2_com-r1_1, ddr2_com-Ie1\ddr1_o ) + m3*cross( ...
    r3_com-r1_1, ddr3_com-Ie1\ddr1_o ) + I2*theDdot2 + I3*theDdot3 ;

% H2_1 = m2*cross(r2_com-r1_1,ddr2_com ) + m3 * cross( r3_com-r1_1, ...
%     ddr3_com ) + I2*theDdot2 + I3*theDdot3 ;

Mo_2 = (cross( (r2_com-r1_1), Fg2 ) + cross( (r3_com-r1_1), Fg3 ) ) ;

r_23_com = ( m2 * r2_com + m3 * r3_com ) / (m2+m3) ;

AMB2 = simplify( H2_1 == Mo_2 + (m2+m3)*cross( r1_1-r_23_com, Ie1\ddr1_o ) ) ;
% AMB2 = simplify( H2_1 == Mo_2 ) ;

% AMB of the third link about joint 2:

% Acceleration of the com of links 1 and 2:

ddr2_o = Ie1\ddr1_o + Ie2\dr2_1 ;

% H3_2 = m3*cross( r3_com-r2_1, ddr3_com-ddr2_o ) + I3*theDdot3 ;
H3_2 = m3 * cross( r3_com-r2_1, ddr3_com ) + I3 * theDdot3 ;

Mo_3 = cross( (r3_com-r2_1), Fg3 ) ;

r12_com = ( m1 * r1_com + m2 * r2_com ) / ( m1 + m2 ) ;

% AMB3 = simplify( H3_2 == Mo_3 + m3*cross( r2_1-r3_com, ddr2_o ) ) ;
AMB3 = simplify( H3_2 == Mo_3 ) ;

% Solve for the EOM:

EOM = solve( [ AMB1(3), AMB2(3), AMB3(3)], [ddthet1, ddthet2, ddthet3 ] ) ;

thetDDot1 = simplify(EOM.ddthet1) ;
thetDDot2 = simplify(EOM.ddthet2) ;
thetDDot3 = simplify(EOM.ddthet3) ;

matlabFunction( thetDDot1, 'File', 'the1DotDot' ) ;
matlabFunction( thetDDot2, 'File', 'the2DotDot' ) ;
matlabFunction( thetDDot3, 'File', 'the3DotDot' ) ;
