function simulateDAE()

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
    
    th10 = -pi/3 ;
    th20 = -pi/3;
    th30 = -pi/3 ;

    dth10 = 0 ;
    dth20 = 0 ;
    dth30 = 0 ;
    
    x10 = p.l1 * cos( th10 ) ;
    x20 = p.l1 * cos( th10 ) + p.l2 * cos( th20 ) ;
    x30 = p.l1 * cos( th10 ) + p.l2 * cos( th20 ) + p.l3 * cos( th30 ) ;
    
    y10 = p.l1 * sin( th10 ) ;
    y20 = p.l1 * sin( th10 ) + p.l2 * sin( th20 ) ;
    y30 = p.l1 * sin( th10 ) + p.l2 * sin( th20 ) + p.l3 * sin( th30 ) ;
    
    dx10 = p.l1 * dth10 * cos( th10 ) ;
    dx20 = dx10 + p.l2 * dth20 * cos( th20 ) ;
    dx30 = dx20 + p.l3 * dth30 * cos( th30 ) ;
    
    dy10 = p.l1 * dth10 * sin( th10 ) ;
    dy20 = dy10 + p.l2 * dth20 * sin( th20 ) ;
    dy30 = dy20 + p.l3 * dth30 * sin( th30 ) ;
    
    z0 = [ x10; x20; x30; y10; y20; y30; th10; th20; th30; ...
        dx10; dx20; dx30; dy10; dy20; dy30; dth10; dth20; dth30 ] ;
    
    % Time span:
    tSpan = 0 : 0.1 : 100 ;
    
    % ODE45 options:
    opts = odeset( 'AbsTol', 1e-10, 'RelTol', 1e-10 ) ;
    
    % ODE45:
    [ time, state ] = ode45( @dyn, tSpan, z0, opts, p ) ;
    
    % Save Thetas:
    filename = 'C:\Users\Jonah\Documents\MAE 4730\MAE 4730\Homework\Final Project\Thetas_To_Compare\thetasDAE' ;
    mat = [state( :, 7 : 9 ), state(:, 16 : end ) ] ;
    save( filename, 'mat' );
    
    % Save Parmeters:
    filename = 'C:\Users\Jonah\Documents\MAE 4730\MAE 4730\Homework\Final Project\Thetas_To_Compare\paramDAE' ;
    mat = p ;
    save( filename, 'mat' ) ;
    
    % Setup simulation:
    
    x1 = state(:,1) ;
    y1 = state(:,4) ;
    
    x2 = state(:,2) ;
    y2 = state(:,5) ;
    
    x3 = state(:,3) ;
    y3 = state(:,6) ;
    
    % Plot trajectories:
    figure( 1 ) ;
    plot( x1, y1 ) ;
    hold on ;
    plot( x2, y2 ) ;
    plot( x3, y3 ) ;
    title( 'DAE Trajectory of Each Link' ) ;
    xlabel( 'X-Position' ) ;
    ylabel( 'Y-Position' ) ;
    hold off ;
    legend( 'Link 1', 'Link 2', 'Link 3' ) ;
    
    % Begin Animation:
    figure( 2 ) ;
    bar1 = plot( [0 x1(1)], [0 y1(1)], 'b', 'LineWidth', 2 ) ;
    hold on ;
    bar2 = plot( [x1(1) x2(1)], [y1(1) y2(1)], 'b', 'LineWidth', 2 ) ;
    bar3 = plot( [x2(1) x3(1)], [y2(1) y3(1)], 'b', 'LineWidth', 2 ) ;
    title( 'DAE Triple Pendulum Animation' ) ;
    xlabel( 'X-Position' ) ;
    ylabel( 'Y-Position' ) ;
    
    axis ([ -5 5 -5 5 ]) ;
    %axis equal ;
    hold off ;
    
    t = 0 ;
    tic ;
    prop = 3 ;
    
    for i = 1 : size( time, 1 )
       
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
        x11 = z( 1 ) ;
        x21 = z( 2 ) ;
        x31 = z( 3 ) ;
        
        y11 = z( 4 ) ;
        y21 = z( 5 ) ;
        y31 = z( 6 ) ;
        
        th1 = z( 7 ) ;
        th2 = z( 8 ) ;
        th3 = z( 9 ) ;
        
        dth1 = z( 16 ) ;
        dth2 = z( 17 ) ;
        dth3 = z( 18 ) ;
        
        % Take derivative term:
        
        A = AMatrix(p.I1,p.I2,p.I3,p.d1,p.d2,p.d3,p.l1,p.l2,p.l3,p.m1,p.m2,p.m3,th1,th2,th3,x11,...
            x21,x31,y11,y21,y31 ) ;
        b = bVector(p.d1,p.d2,p.d3,dth1,dth2,dth3,p.g,p.l1,p.l2,p.m1,p.m2,p.m3,th1,...
            th2,th3) ;
        
        ddz = A\ b ;
        
        dz = [ z(10:18); ddz(1:9) ] ;
        
    end

end