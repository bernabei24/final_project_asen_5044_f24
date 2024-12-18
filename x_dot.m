function dxdt = x_dotODE45(t,x, u_func,L)
% state_dependent_system computes the derivative of x for a state-dependent system
    % Inputs:
    %   t - current time (required by ode45 but not used here)
    %   x - current state vector (6x1)
    %   u - constant control input vector (4x1)
    % Output:
    %   dxdt - derivative of state vector (6x1)
    u = u_func(t,x);
    % Define the nonlinear dynamics for each state
    
    dx1 = u(1) * cos(x(3)); % Easting of ground
    dx2 = u(1) * sin(x(3));  % Example for x2 (Northing of ground)
    dx3 = (u(1)/L) * tan(u(2));  % Example for x3 (Heading of ground)
    dx4 = u(3) * cos(x(6));   % Example for x4 (Easting of air)
    dx5 = u(3) * sin(x(6));   % Example for x5 (Northing of air)
    dx6 = u(4);               % Example for x6 (Heading of air)

    % Combine all the derivatives into one vector
    dxdt = [dx1; dx2; dx3; dx4; dx5; dx6];
end