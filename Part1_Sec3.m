clear;
clc;
close all;

%initial conditions

L = 0.5; %m 
Eg_init = 10; %m
Ng_init = 0; %m
thetag_init = pi/2; %rad
vg_init = 2; %m/s
phi_init = -pi/18; %rad
Ea_init = -60; %m
Na_init = 0; %m
thetaa_init = -pi/2; %rad
va_init = 12; %m/s
omegaa_init = pi/25; %rad/s
dt = 0.1; %sec

%constructing initial vectors
perturb_x0 = [0; 1; 0; 0; 0; 0.1];
x_init = [Eg_init Ng_init thetag_init Ea_init Na_init thetaa_init]';
u_init = [vg_init phi_init va_init omegaa_init]';

%Define A(x,u)
A_func = @(x,u) [0 0 -u(1,1)*sin(x(3,1)) 0 0 0; %not necessary for ODE45, but may be helpful for linearization
             0 0  u(1,1)*cos(x(3,1)) 0 0 0;
             0 0  0                  0 0 0;
             0 0 0 0 0  u(3,1)*sin(x(6,1));
             0 0 0 0 0 -u(3,1)*cos(x(6,1));
             0 0 0 0 0 0];

%Define B(x,u)
B_func = @(x,u) [cos(x(3,1)) 0 0 0; %not necessary for ODE45, but may be helpful for linearization
             sin(x(3,1)) 0 0 0;
             tan(u(2,1))/L (u(1,1)/L)*( sec(u(2,1))^2 ) 0 0;
             0 0 cos(x(6,1)) 0;
             0 0 sin(x(6,1)) 0;
             0 0 0 1];
%Define input
u_func = @(t, x) u_init; % Constant control input

%Timespan for simulation
t_span = 0:0.1:100; % Time from 0 to 100 seconds in 0.1-second steps

% Define the nonlinear dynamics
dynamics = @(t, x) x_dotODE45(t, x, u_func,L);

% Solve using ode45
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PERTURBATION ADDED
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[t, x] = ode45(dynamics, t_span, x_init + perturb_x0); %perturbation added

% Wrap the angles theta_g (x(3,:)) and theta_a (x(6,:)) to [-pi, pi]
x(:, 3) = mod(x(:, 3) + pi, 2*pi) - pi;  % Wrap theta_g (ground heading)
x(:, 6) = mod(x(:, 6) + pi, 2*pi) - pi;  % Wrap theta_a (air heading)

% Plot results
figure;
for i = 1:size(x, 2)
    subplot(size(x, 2), 1, i);
    plot(t, x(:, i), 'LineWidth', 1.5);
    grid on;
    
    % Adjust titles and labels based on the state variables
    if i == 1
        title('$\xi$ (Easting of ground)', 'Interpreter', 'latex');
        ylabel('$\xi$ (m)', 'Interpreter', 'latex');
    elseif i == 2
        title('$\eta$ (Northing of ground)', 'Interpreter', 'latex');
        ylabel('$\eta$ (m)', 'Interpreter', 'latex');
    elseif i == 3
        title('$\theta$ (Heading of ground)', 'Interpreter', 'latex');
        ylabel('$\theta$ (rad)', 'Interpreter', 'latex');
    elseif i == 4
        title('$\xi$ (Easting of air)', 'Interpreter', 'latex');
        ylabel('$\xi$ (m)', 'Interpreter', 'latex');
    elseif i == 5
        title('$\eta$ (Northing of air)', 'Interpreter', 'latex');
        ylabel('$\eta$ (m)', 'Interpreter', 'latex');
    elseif i == 6
        title('$\theta$ (Heading of air)', 'Interpreter', 'latex');
        ylabel('$\theta$ (rad)', 'Interpreter', 'latex');
    end
    xlabel('Time (s)', 'Interpreter', 'latex');
end
sgtitle('State Trajectories over Time With Perturbation', 'Interpreter', 'latex');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulate using LTV DT Model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% solve nominal trajectory with no perturbations
[t, x] = ode45(dynamics, t_span, x_init);

% initial conditions for perturbation state vector
x_sim = perturb_x0';

% simulate system behavior at each time step k
for k = 1:1:1000
    % Use CT jacobians to find A and B at t = t_k and nom[k]
    A_nom_k = get_Abar(u_init, x(k, :)');
    B_nom_k = get_Bbar(u_init, x(k, :)', L);

    % Eulerized estimate of DT Jacobians
    F_k = dt * A_nom_k;
    G_k = dt * B_nom_k;

    % Calculate perturbation state at t = t_k
    x_k = F_k * x_sim(k, :)' + G_k * u_init;
    x_sim = [x_sim; x_k'];
end

% find total state by combining nominal trajectory with perturbation states
x_total = x + x_sim;

% Wrap the angles theta_g (x_total(3,:)) and theta_a (x_total(6,:)) to [-pi, pi]
x_total(:, 3) = mod(x_total(:, 3) + pi, 2*pi) - pi;  % Wrap theta_g (ground heading)
x_total(:, 6) = mod(x_total(:, 6) + pi, 2*pi) - pi;  % Wrap theta_a (air heading)

% Plot total state results
figure;
for i = 1:size(x_total, 2)
    subplot(size(x_total, 2), 1, i);
    plot(t, x_total(:, i), 'LineWidth', 1.5);
    grid on;
    
    % Adjust titles and labels based on the state variables
    if i == 1
        title('$\xi$ (Easting of ground)', 'Interpreter', 'latex');
        ylabel('$\xi$ (m)', 'Interpreter', 'latex');
    elseif i == 2
        title('$\eta$ (Northing of ground)', 'Interpreter', 'latex');
        ylabel('$\eta$ (m)', 'Interpreter', 'latex');
    elseif i == 3
        title('$\theta$ (Heading of ground)', 'Interpreter', 'latex');
        ylabel('$\theta$ (rad)', 'Interpreter', 'latex');
    elseif i == 4
        title('$\xi$ (Easting of air)', 'Interpreter', 'latex');
        ylabel('$\xi$ (m)', 'Interpreter', 'latex');
    elseif i == 5
        title('$\eta$ (Northing of air)', 'Interpreter', 'latex');
        ylabel('$\eta$ (m)', 'Interpreter', 'latex');
    elseif i == 6
        title('$\theta$ (Heading of air)', 'Interpreter', 'latex');
        ylabel('$\theta$ (rad)', 'Interpreter', 'latex');
    end
    xlabel('Time (s)', 'Interpreter', 'latex');
end
sgtitle('States vs. Time, Linearized Approximate Dynamics Simulation', 'Interpreter', 'latex');

% plot perturbation state results
% Wrap the angles theta_g (x_total(3,:)) and theta_a (x_total(6,:)) to [-pi, pi]
x_sim(:, 3) = mod(x_sim(:, 3) + pi, 2*pi) - pi;  % Wrap theta_g (ground heading)
x_sim(:, 6) = mod(x_sim(:, 6) + pi, 2*pi) - pi;  % Wrap theta_a (air heading)

% Plot results
figure;
for i = 1:size(x_total, 2)
    subplot(size(x_sim, 2), 1, i);
    plot(t, x_sim(:, i), 'LineWidth', 1.5);
    grid on;
    
    % Adjust titles and labels based on the state variables
    if i == 1
        title('$\xi$ (Easting of ground)', 'Interpreter', 'latex');
        ylabel('$\xi$ (m)', 'Interpreter', 'latex');
    elseif i == 2
        title('$\eta$ (Northing of ground)', 'Interpreter', 'latex');
        ylabel('$\eta$ (m)', 'Interpreter', 'latex');
    elseif i == 3
        title('$\theta$ (Heading of ground)', 'Interpreter', 'latex');
        ylabel('$\theta$ (rad)', 'Interpreter', 'latex');
    elseif i == 4
        title('$\xi$ (Easting of air)', 'Interpreter', 'latex');
        ylabel('$\xi$ (m)', 'Interpreter', 'latex');
    elseif i == 5
        title('$\eta$ (Northing of air)', 'Interpreter', 'latex');
        ylabel('$\eta$ (m)', 'Interpreter', 'latex');
    elseif i == 6
        title('$\theta$ (Heading of air)', 'Interpreter', 'latex');
        ylabel('$\theta$ (rad)', 'Interpreter', 'latex');
    end
    xlabel('Time (s)', 'Interpreter', 'latex');
end
sgtitle('Linearized Approx Perturbations vs. Time', 'Interpreter', 'latex');





