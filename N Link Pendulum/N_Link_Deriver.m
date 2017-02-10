function N_Link_Deriver(N)

%N = 4 ;

% unit vectors ;
i = [ 1; 0; 0 ] ;
j = [ 0; 1; 0 ] ;
k = [ 0; 0; 1 ] ;

% Constants needed:
syms g real ;

% Variables needed:
th      = sym( 'th', [ N 1 ], 'real' ) ;
dth     = sym( 'dth', [ N 1 ], 'real' ) ;
ddth    = sym( 'ddth', [ N 1 ], 'real' ) ;
l       = sym( 'l', [ N 1 ], 'real' ) ;
m       = sym( 'm', [ N 1 ], 'real' ) ;
I       = sym( 'I', [ N 1 ], 'real' ) ;
AMB     = sym( zeros( N, 1) ) ;
H       = sym( zeros( N, 1) ) ;

Fg      = -m * g * j' ; 
thDot   = dth * k' ; 
thDdot  = ddth * k' ;

IeN = [ cos(th), sin(th), 0 ; ...
       -sin(th), cos(th), 0 ; ...
       0, 0, 1 ] ;
   
r      = l * i' ;
dr     = cross( thDot, r ) ;
ddr    = sym( zeros(N,3) ) ;
r_o    = sym( zeros(N,3) ) ;
dr_o   = sym( zeros(N,3) ) ;
ddr_o  = sym( zeros(N,3) ) ; 

r_com   = sym( zeros( N, 3) ) ;
a_com   = sym( zeros( N, 3) ) ;

ddr(1,:) = simplify(jacobian( dr(1,:), dth(1) ) * ddth(1) + ...
    cross( thDot(1,:), dr(1,:), 2 )' ) ;

r_o(1,:)    = simplify( IeN([1 (N+1) end], : )\ r( 1, : )' ) ;
r_com(1,:)  = simplify( r_o(1,:)/ 2 ) ;
dr_o(1,:)   = simplify( IeN( [1 (N+1) end], : )\ dr( 1, : )' ) ;
ddr_o(1,:)  = simplify( IeN( [1 (N+1) end], : )\ ddr( 1, : )' ) ;
a_com(1,:)  = simplify( ddr_o(1,:)/ 2 ) ;

for i = 2 : N 

    ddr(i,:)    = jacobian( dr(i,:), dth(i) ) * ddth(i) + cross( thDot(i,:), dr(i,:) )' ;
    r_o(i,:)    = simplify( r_o(i-1,:) + (IeN( [i (N+i) end], : ) ...
        \ r( i, : )')' ) ;
    dr_o(i,:)   = simplify( dr_o(i-1,:) + (IeN( [i (N+i) end], : )...
        \ dr(i, :)')' ) ;
    ddr_o(i,:)  = simplify( ddr_o(i-1,:) + (IeN( [i (N+i) end], : )...
        \ ddr( i, : )')' ) ;
    
    r_com(i,:)  = simplify( r_o(i-1,:) + (IeN([i (N+i) end],:)\ r(i,:)'/2)' ) ;
    a_com(i,:)  = simplify( ddr_o(i-1,:) + (IeN([i (N+i) end],:)\ddr(i,:)'/ 2)' ) ;

end

i = 2 ;

H(1) = dot( simplify( sum( m' * cross( r_com, a_com, 2 ) , 1 ) + I' * thDdot ), k ) ;
M = dot( sum( simplify( cross( r_com,  Fg, 2 ) ), 1 ), k ) ;
AMB(1) = simplify( H(1) == M ) ;

while i < N
   
    %H_i = simplify( dot( sum( cross( r_com( i : end, : ), a_com( i : end, : ) ), 1 ), k) );
    
    H_i = ones(N-i+1,1)*r_o(i-1,:) ;
    
    H(i) = dot( simplify( sum( m( i : end )' * cross( r_com( i : end, : )-H_i, ...
        a_com( i : end, : ), 2 ), 1 ) + I( i : end )' * thDdot( i : end, : ) ), k ) ; 
    
    r_eff = r_com( i : end, : ) ;
    F_i = Fg( i : end, : ) ;
    
    rxf = dot( sum( cross( -(H_i-r_eff), F_i, 2 ), 1 ), k ) ;
    
    %M   = dot( simplifY( sum(cross( r_com( i : end, : ) - ones(N-i+1,1)*r_o(i-1,:), ...
     %   Fg( i : end, : ) ), 1 ) ), k ) ;
    %M = dot( simplify( (H_i-r_eff)' * Fg(i : end, : ) ), k ) ;
    
    AMB( i ) = simplify( H(i) == rxf ) ;
    
    i = i + 1 ;
    
end

H(end) = dot( cross( r_com(end,:)-r_o(end-1,:), a_com(end,:) ) + ...
    I(end)'*thDdot(end,:), k ) ;
AMB( end ) = simplify( H(end) == dot( sum(cross(r_com(end,:)-r_o(end-1,:),...
    Fg(end,:),2 ), 1 ) ,k) ) ;

% EOM = solve( AMB, ddth ) ;
% fields = fieldnames( EOM ) ;
% 
% eqns = sym( zeros( numel(fields), 1 ) ) ;
% 
% for i = 1 : numel(fields)
%    
%     eqns( i ) = EOM.(fields{i}) ;
%     
% end

[ A, b ] = equationsToMatrix( AMB, ddth ) ;
matlabFunction( A, 'File', 'A_matrix', 'Vars', {I', l', m', th'} ) ;
matlabFunction( b, 'File', 'b_matrix', 'Vars', {dth', g, l', m', th'}  ) ;

end