function N_Link_Simulator() 
    
    N_prev = open( 'N.mat' ) ;

    p.N = 9 ;
    N = p.N ;
    
    save( 'N.mat' ) ;

    if ( N_prev.N ~= p.N ) 
       
        N_Link_Deriver( p.N ) ;
        
    end
    
    % Constant parameters
    p.g = 1 ;
    p.l = ones( p.N, 1 ) ;
    p.m = ones( p.N, 1 ) ;
    p.I = ones( p.N, 1 ) ;
    
    % Initial conditions:
    the1_0 = -pi/3 ;
    the2_0 = -pi/3 ;
    the3_0 = -pi/3 ;
    
    th0 = ones( p.N, 1 );
    dth0 = zeros( p.N, 1 ) ;
    
    for i = 1 : p.N 
       
        th0(i) = pi/3 + pi/100 * (i/p.N) ;
        p.I(i) = 1/12 * p.m(i) * p.l(i)^2 ;
        
    end
    
    dthe1_0 = 0 ;
    dthe2_0 = 0 ;
    dthe3_0 = 0 ;
    
    %z0 = [ the1_0; the2_0; the3_0; dthe1_0; dthe2_0; dthe3_0 ] ;
    z0 = [ th0; dth0 ];
    
    % Time span:
    tSpan = [ 0 100 ] ;
    
    % ODE45 options:
    opts = odeset( 'AbsTol', 1e-10, 'RelTol', 1e-10 ) ;
    
    % ODE45:
    [ time, state ] = ode45( @dyn, tSpan, z0, opts, p ) ;
    
    % Save Theta Values:
%     filename = 'C:\Users\Jonah\Documents\MAE 4730\MAE 4730\Homework\Final Project\Thetas_To_Compare\thetasAMB' ;
%     mat = state( :, 1 : 3 ) ;
%     save( filename, 'mat' );
    
    % Setup simulation:
    
%     x1 = p.l1 * cos( state(:,1) ) ;
%     y1 = p.l1 * sin( state(:,1) ) ;
%     
%     x2 = p.l2 * cos( state(:,2) ) + x1 ;
%     y2 = p.l2 * sin( state(:,2) ) + y1 ;
%     
%     x3 = p.l3 * cos( state(:,3) ) + x2 ;
%     y3 = p.l3 * sin( state(:,3) ) + y2 ;
    
    x = ones( size(state,1), p.N ) ;
    y = ones( size(state,1), p.N ) ;
    
    x(:,1) = p.l(1) * cos( state(:,1) ) ;
    y(:,1) = p.l(1) * sin( state(:,1) ) ;
    
    for i = 2 : p.N
       
        x(:,i) = p.l(i) * cos( state(:,i) ) + x(:,i-1) ;
        y(:,i) = p.l(i) * sin( state(:,i) ) + y(:,i-1) ;
        
    end
    
%     figure( 1 ) ;
%     plot( x2, y2 ) ;
    
    figure( 1 ) ;
%     bar1 = plot( [0 x1(1)], [0 y1(1)], 'b', 'LineWidth', 2 ) ;
%     hold on ;
%     bar2 = plot( [x1(1) x2(1)], [y1(1) y2(1)], 'b', 'LineWidth', 2 ) ;
%     bar3 = plot( [x2(1) x3(1)], [y2(1) y3(1)], 'b', 'LineWidth', 2 ) ;
    
    bar = plot( [ 0 x(1,:)], [ 0 y(1,:)], 'b', 'LineWidth', 3 ) ;
    hold on ;
    
    axis ([ -10 10 -10 10 ]) ;
    %axis equal ;
    hold off ;
    
    t = 0 ;
    tic ;
    prop = 5 ;
    
    for i = 1 : size( time, 1 )
       
%         x1_1 = interp1( time, x1, prop*t ) ;
%         y1_1 = interp1( time, y1, prop*t ) ;
%         
%         x2_1 = interp1( time, x2, prop*t ) ;
%         y2_1 = interp1( time, y2, prop*t ) ;
%         
%         x3_1 = interp1( time, x3, prop*t ) ;
%         y3_1 = interp1( time, y3, prop*t ) ;
        
%         set( bar1, 'XData', [0 x1_1], 'YData', [0 y1_1] ) ;
%         set( bar2, 'XData', [x1_1 x2_1], 'YData', [y1_1 y2_1] ) ;
%         set( bar3, 'XData', [x2_1 x3_1], 'YData', [y2_1 y3_1] ) ;

        x_1 = interp1( time, x, prop*t );
        y_1 = interp1( time, y, prop*t );
        
        set( bar, 'XData', [0 x_1], 'YData', [ 0 y_1] ) ;
    
        drawnow ;
        
        t = toc ;
        
    end
    
    function dz = dyn( ~, z, p ) 
        
        A = A_matrix( p.I', p.l', p.m', z(1:p.N)' ) ;
        b = b_matrix( z((p.N+1):end)', p.g, p.l', p.m', z(1:p.N)' ) ;
        
        % Take the derivative: 
        
        ddth = A\ b ;
        
        dz = [ z((p.N+1):end); ddth ] ;
        
    end

end