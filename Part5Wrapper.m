% Extended Kalman Filter
clear;
clc;

load("cooplocalization_finalproj_KFdata.mat")

%GENERATE RANDOM PROCESS NOISE VECTORS
% Check if all eigenvalues of Q are positive
eigenvaluesQ = eig(Qtrue);
rng(1); % for reproduceability
if all(eigenvaluesQ > 0)
    disp('The matrix Q is positive definite.');

    % Perform the Cholesky decomposition
    LowerQ = chol(Qtrue, 'lower');

    % Generate a matrix of standard normal random variables
    
    ZQ = randn(6, 1001);

    % Multiply by the Cholesky factor to get the random vectors with covariance Q
    w_k = LowerQ * ZQ;

else
    disp('The matrix Q is not positive definite.');
end

%GENERATE RANDOM MEASUREMENT NOISE VECTORS
eigenvaluesR = eig(Rtrue);
if all(eigenvaluesR > 0)
    disp('The matrix R is positive definite.');

    % Perform the Cholesky decomposition
    LowerR = chol(Rtrue, 'lower');

    % Generate a matrix of standard normal random variables
    
    ZR = randn(5, 1001);

    % Multiply by the Cholesky factor to get the random vectors with covariance Q
    v_k = LowerR * ZR;

else
    disp('The matrix R is not positive definite.');
end

u_func = @(t, x) u_init; % Constant control input
w_func = @(t) w_k(:, floor(t / 0.1) + 1); % 

% Define the nonlinear dynamics
dynamics_noise = @(t, x) x_dotODE45noise(t, x, u_func, w_func, L); %for truth model
dynamics_nominal = @(t, x) x_dotODE45(t, x, u_func,L); %for nominal trajectory

% step 1: define initial conditions

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

perturb_x0 = [0; 1; 0; 0; 0; 0.1];
x_init = [Eg_init Ng_init thetag_init Ea_init Na_init thetaa_init]';
u_init = [vg_init phi_init va_init omegaa_init]';

% Define the nonlinear dynamics
u_func = @(t, x) u_init; % Constant control input
w_func = @(t) w_k(:, floor(t / 0.1) + 1); % 
dynamics_noise = @(t, x) x_dotODE45noise(t, x, u_func, w_func, L); %for truth model
dynamics_nominal = @(t, x) x_dotODE45(t, x, u_func,L); %for nominal trajectory

% solve nominal trajectory with no perturbations
[t, x_true] = ode45(dynamics_noise, tvec, x_init);  %perturbation not added
[t, x_nominal] = ode45(dynamics_nominal, tvec, x_init); %perturbation not added
y_true = get_Y(x_true);
y_nominal = get_Y(x_nominal);
% define ydata for testing
% ydata = ydata;    % default uses given ydata from canvas
% ydata = y_true';

% x_init = [Eg_init Ng_init thetag_init Ea_init Na_init thetaa_init]' + perturb_x0;

P_init = 1000 * diag([1, 1, 1, 1, 1, 1]);

%pretty sure this is what gamma looks like although not positive
Gamma = eye(6);
Omega = dt * Gamma;

% set up CT nonlinear dynamics functions
f = @(t, x) x_dotODE45(t, x, u_func, L); % for computing trajectory online

% step 2: k = 0
x_min = x_init;
x_plus = x_init;
P_min = P_init;
P_plus = P_init;

for k = 1:1000
    % step 3: time update/prediction step for time k + 1
    x_plus_k = x_plus(:, k);
    [t, x_min_kplus1] = ode45(f, tvec(k:k+1), x_plus_k);
    x_min(:, k + 1) = x_min_kplus1(end, :)';

    A = get_Abar(u_init, x_plus(:, k));
    F_k = eye(6) + dt * A;
    P_min(:, :, k + 1) = F_k * P_plus(:, :, k) * F_k' + Omega * Qtrue * Omega';

    % step 4: measurement update/correction step for k + 1
    y_min = get_Y_k(x_min(:, k + 1));
    H = get_Cbar(x_min(:, k + 1));
    e_y = ydata(:, k + 1) - y_min;
    K = P_min(:, :, k+1) * H' * inv(H * P_min(:, :, k+1) * H' + Rtrue);

    x_plus(:, k + 1) = x_min(:, k + 1) + K * e_y;
    P_plus(:, :, k + 1) = (eye(6) - K * H) * P_min(:, :, k+1);
end

x_plus = x_plus';
t = dt*(0:1000);

% plot estimated state results
% Wrap the angles theta_g (x_total(3,:)) and theta_a (x_total(6,:)) to [-pi, pi]
x_plus(:, 3) = mod(x_plus(:, 3) + pi, 2*pi) - pi;  % Wrap theta_g (ground heading)
x_plus(:, 6) = mod(x_plus(:, 6) + pi, 2*pi) - pi;  % Wrap theta_a (air heading)

EKF_total_state_graphs(x_nominal, x_plus', P_plus, tvec)

% % Plot results
% figure;
% for i = 1:size(x_plus, 2)
%     subplot(size(x_plus, 2), 1, i);
%     plot(t, x_plus(:, i), 'LineWidth', 1.5);
%     grid on;
% 
%     % Adjust titles and labels based on the state variables
%     if i == 1
%         title('$\xi$ (Easting of ground)', 'Interpreter', 'latex');
%         ylabel('$\xi$ (m)', 'Interpreter', 'latex');
%     elseif i == 2
%         title('$\eta$ (Northing of ground)', 'Interpreter', 'latex');
%         ylabel('$\eta$ (m)', 'Interpreter', 'latex');
%     elseif i == 3
%         title('$\theta$ (Heading of ground)', 'Interpreter', 'latex');
%         ylabel('$\theta$ (rad)', 'Interpreter', 'latex');
%     elseif i == 4
%         title('$\xi$ (Easting of air)', 'Interpreter', 'latex');
%         ylabel('$\xi$ (m)', 'Interpreter', 'latex');
%     elseif i == 5
%         title('$\eta$ (Northing of air)', 'Interpreter', 'latex');
%         ylabel('$\eta$ (m)', 'Interpreter', 'latex');
%     elseif i == 6
%         title('$\theta$ (Heading of air)', 'Interpreter', 'latex');
%         ylabel('$\theta$ (rad)', 'Interpreter', 'latex');
%     end
%     xlabel('Time (s)', 'Interpreter', 'latex');
% end
% sgtitle('EKF Estimated States vs. Time', 'Interpreter', 'latex');