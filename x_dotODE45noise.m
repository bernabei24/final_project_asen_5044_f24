function dxdt = x_dotODE45noise(t,x, u_func,w_func, L)
% state_dependent_system computes the derivative of x for a state-dependent system
    % Inputs:
    %   t - current time (required by ode45 but not used here)
    %   x - current state vector (6x1)
    %   u - constant control input vector (4x1)
    %   w_func — gets predetermined process noise at time t
    % Output:
    %   dxdt - derivative of state vector (6x1)
    u = u_func(t,x);
    w = w_func(t);
    % Define the nonlinear dynamics for each state
    
    dx1 = u(1) * cos(x(3)) + w(1); % Easting of ground
    dx2 = u(1) * sin(x(3)) + w(2);  % Example for x2 (Northing of ground)
    dx3 = (u(1)/L) * tan(u(2)) + w(3);  % Example for x3 (Heading of ground)
    dx4 = u(3) * cos(x(6)) + w(4);   % Example for x4 (Easting of air)
    dx5 = u(3) * sin(x(6)) + w(5);   % Example for x5 (Northing of air)
    dx6 = u(4) + w(6);               % Example for x6 (Heading of air)

    % Combine all the derivatives into one vector
    dxdt = [dx1; dx2; dx3; dx4; dx5; dx6];
end