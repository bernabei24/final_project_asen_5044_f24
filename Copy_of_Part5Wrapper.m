% Extended Kalman Filter

%%%%%%%%%%%%%%%%%%%%%%%%
% please be careful messing around with this copy. This is the version of
% the EKF that exhibits the best performance so far.
%%%%%%%%%%%%%%%%%%%%%%%%

clear;
clc;
close all;
load("cooplocalization_finalproj_KFdata.mat")

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

%constructing initial vectors
perturb_x0 = [0; 1; 0; 0; 0; 0.1];
x_init = [Eg_init Ng_init thetag_init Ea_init Na_init thetaa_init]' + perturb_x0;
u_init = [vg_init phi_init va_init omegaa_init]';

%Define input
u_func = @(t, x) u_init; % Constant control input
w_func = @(t) w_k(:, floor(t / 0.1) + 1);

% Define the nonlinear dynamics
dynamics_noise = @(t, x) x_dotODE45noise(t, x, u_func, w_func, L); %for truth model
dynamics_nominal = @(t, x) x_dotODE45(t, x, u_func,L); %for nominal trajectory

% solve nominal trajectory with no perturbations
[t, x_true] = ode45(dynamics_noise, tvec, x_init);
[t, x_nominal] = ode45(dynamics_nominal, tvec, x_init);

% Wrap the angles theta_g (x(3,:)) and theta_a (x(6,:)) to [-pi, pi]
x_true(:, 3) = mod(x_true(:, 3) + pi, 2*pi) - pi;  % Wrap theta_g (ground heading)
x_true(:, 6) = mod(x_true(:, 6) + pi, 2*pi) - pi;  % Wrap theta_a (air heading)
x_nominal(:, 3) = mod(x_nominal(:, 3) + pi, 2*pi) - pi;  % Wrap theta_g (ground heading)
x_nominal(:, 6) = mod(x_nominal(:, 6) + pi, 2*pi) - pi;  % Wrap theta_a (air heading)

% find total measurement vector
y_true = get_Y_noise(x_true,v_k');
%y_true(:, 1) = mod(y_true(:, 1) + pi, 2*pi) - pi;  % Wrap gamma_ag
%y_true(:, 3) = mod(y_true(:, 3) + pi, 2*pi) - pi;  % Wrap gamme_ga
y_nom = get_Y(x_nominal);

% define ydata for testing
% ydata = ydata;    % default uses given ydata from canvas
ydata = y_true';

%%%%%%%
% INITIALIZING EKF LOOP
P_init = 0.1 * diag([1, 1, 1, 1, 1, 1]);
%%%%%%%
% gonna try and tune Q in here
% 0.6*Qtrue seems to perform the best
%Qtrue = 0.6*Qtrue;
%%%%%%%

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
    %P_min(:, :, k + 1) = F_k * P_plus(:, :, k) * F_k' + Qtrue;


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
x_true(:, 3) = mod(x_true(:, 3) + pi, 2*pi) - pi;  % Wrap theta_g (ground heading)
x_true(:, 6) = mod(x_true(:, 6) + pi, 2*pi) - pi;  % Wrap theta_a (air heading)

% plots with 2sigma error bounds
%EKF_total_state_graphs(x_true, x_plus', P_plus, tvec)
plot_states_and_errorsEKF(x_plus,x_true,P_plus,tvec)

