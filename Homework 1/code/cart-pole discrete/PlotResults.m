function MovieFile = PlotResults(sol)
    time = sol.time;
    state = sol.state;
    control = sol.control;
    cost = sol.cost;
    target = sol.target;
    K = sol.gain;
    l = sol.l;

    % Plot Results %
    figure('Position',[300 100 624 564]);

    % Plot States: x [m], v [m/sec], theta [deg], thetadot [deg/sec]
    plot(state(1, :), time, 'k');   % x [m]
    title('Cart x-position over time')
    xlabel('x-position [m]')
    ylabel('Time [s]')
    
    figure;
    plot(time, state(2, :), 'b');   % v [m/sec]
    title('Cart velocity, v [m/sec]')
    xlabel('x-velocity [m/sec]')
    ylabel('Time [s]')
    
    figure;
    plot(time, state(3, :), 'r');   % theta [deg]
    title('Cart rod angle')
    ylabel('Rod angle, theta [deg]')
    xlabel('Time [s]')
    
    figure;
    plot(time, state(4, :), 'g');   % thetadot [deg/sec]
    title('Cart rod angular velocity over time')
    ylabel('Rod angular velocity, thetadot [deg/sec]')
    xlabel('Time [s]')
    
    % Plot cost vs number of iterations
    figure;
    plot(1:200, cost);
    title('Cartpole swing-up task cost over time')
    xlabel('Time [s]')
    ylabel('Cost')
    
    % Plot controller gains
    figure('Position',[600 100 524 564]);
    plot(l);
    title('l-gain')
    xlabel('Timestep')
    ylabel('l-gain')

    %  animation
    ...

end







      
