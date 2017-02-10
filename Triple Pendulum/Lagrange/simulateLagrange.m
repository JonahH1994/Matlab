function simulateLagrange() 
    
    clear ;
    close ;
    
    % Constant parameters
    p.g = 1 ;
    p.l1 = 1 ;
    p.l2 = 1 ;
    p.l3 = 1 ;
    p.m1 = 1 ;
    p.m2 = 1 ;
    p.m3 = 1 ;
    
    p.d1 = 0.5 ;
    p.d2 = 0.5 ;
    p.d3 = 0.5 ;
    
    p.I1 = 1/12 * p.m1 * p.l1^2 + p.m1 * (p.d1-p.l1/2).^2;
    p.I2 = 1/12 * p.m2 * p.l2^2 + p.m2 * (p.d2-p.l2/2).^2;
    p.I3 = 1/12 * p.m3 * p.l3^2 + p.m3 * (p.d3-p.l3/2).^2;
    
    % Initial conditions:
    the1_0 = -pi/3 ;
    the2_0 = -pi/3 ;
    the3_0 = -pi/3 ;
    
    dthe1_0 = 0 ;
    dthe2_0 = 0 ;
    dthe3_0 = 0 ;
    
    z0 = [ the1_0; the2_0; the3_0; dthe1_0; dthe2_0; dthe3_0 ] ;
    
    % Time span:
    tSpan = 0 : 0.1 : 100 ;
    
    % ODE45 options:
    opts = odeset( 'AbsTol', 1e-10, 'RelTol', 1e-10 ) ;
    
    % ODE45:
    [ time, state ] = ode45( @dyn, tSpan, z0, opts, p ) ;
    
    % Save thetas: 
    
    filename = 'C:\Users\Jonah\Documents\MAE 4730\MAE 4730\Homework\Final Project\Thetas_To_Compare\thetasLagrange' ;
    mat = state( :, 1 : end ) ;
    save( filename, 'mat' );
    
    % Save Parmaeters:
    filename = 'C:\Users\Jonah\Documents\MAE 4730\MAE 4730\Homework\Final Project\Thetas_To_Compare\paramLagrange' ;
    mat = p ;
    save( filename, 'mat' ) ;
    
    % Setup simulation:
    
    x1 = p.l1 * cos( state(:,1) ) ;
    y1 = p.l1 * sin( state(:,1) ) ;
    
    x2 = p.l2 * cos( state(:,2) ) + x1 ;
    y2 = p.l2 * sin( state(:,2) ) + y1 ;
    
    x3 = p.l3 * cos( state(:,3) ) + x2 ;
    y3 = p.l3 * sin( state(:,3) ) + y2 ;
    
    % Plot trajectories:
    figure( 1 ) ;
    plot( x1, y1 ) ;
    hold on ;
    plot( x2, y2 ) ;
    plot( x3, y3 ) ;
    title( 'Lagrange Trajectory of Each Link' ) ;
    xlabel( 'X-Position' ) ;
    ylabel( 'Y-Position' ) ;
    legend( 'Link 1', 'Link 2', 'Link 3' ) ;
    hold off ;
    
    % Begin Animaiton:
    figure( 2 ) ;
    bar1 = plot( [0 x1(1)], [0 y1(1)], 'b', 'LineWidth', 2 ) ;
    hold on ;
    bar2 = plot( [x1(1) x2(1)], [y1(1) y2(1)], 'b', 'LineWidth', 2 ) ;
    bar3 = plot( [x2(1) x3(1)], [y2(1) y3(1)], 'b', 'LineWidth', 2 ) ;
    title( 'Lagrange Triple Pendulum Animation' ) ;
    xlabel( 'X-Position' ) ;
    ylabel( 'Y-Position' ) ;
    hold off ;
    
    axis ([ -5 5 -5 5 ]) ;
    %axis equal ;
    hold off ;
    
    t = 0 ;
    tic ;
    prop = 3 ;
    
    while ( t < time( end ) )
       
        x1_1 = interp1( time, x1, prop*t ) ;
        y1_1 = interp1( time, y1, prop*t ) ;
        
        x2_1 = interp1( time, x2, prop*t ) ;
        y2_1 = interp1( time, y2, prop*t ) ;
        
        x3_1 = interp1( time, x3, prop*t ) ;
        y3_1 = interp1( time, y3, prop*t ) ;
        
        set( bar1, 'XData', [0 x1_1], 'YData', [0 y1_1] ) ;
        set( bar2, 'XData', [x1_1 x2_1], 'YData', [y1_1 y2_1] ) ;
        set( bar3, 'XData', [x2_1 x3_1], 'YData', [y2_1 y3_1] ) ;
        drawnow ;
        
        t = toc ;
        
    end
    
    function dz = dyn( ~, z, p ) 
       
        % Unpack variables:
        th1 = z( 1 ) ;
        th2 = z( 2 ) ;
        th3 = z( 3 ) ;
        
        dth1 = z( 4 ) ;
        dth2 = z( 5 ) ;
        dth3 = z( 6 ) ;
        
        % Take derivative term:
        ddth1 = the1DotDot(p.I1,p.I2,p.I3,p.d1,p.d2,p.d3,dth1,dth2,dth3,p.g,p.l1,p.l2,p.m1,p.m2, ...
            p.m3, th1, th2, th3 ) ;
        ddth2 = the2DotDot(p.I1,p.I2,p.I3,p.d1,p.d2,p.d3,dth1,dth2,dth3,p.g,p.l1,p.l2,p.m1,p.m2, ...
            p.m3, th1, th2, th3 ) ;
        ddth3 = the3DotDot(p.I1,p.I2,p.I3,p.d1,p.d2,p.d3,dth1,dth2,dth3,p.g,p.l1,p.l2,p.m1,p.m2, ...
            p.m3, th1, th2, th3 ) ;
        
        dz = [ dth1; dth2; dth3; ddth1(1); ddth2(1); ddth3(1) ] ;
        
    end

end