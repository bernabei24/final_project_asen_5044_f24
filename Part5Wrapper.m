% Extended Kalman Filter
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
rng(10); % for reproduceability
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
%perturb_x0 = [0; 1; 0; 0; 0; 0.1];
%nick's test

perturb_x0 = zeros(6,1);
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

% wrapping x_true causes errors in the EKF
x_plot = x_true;

% Wrap the angles theta_g (x(3,:)) and theta_a (x(6,:)) to [-pi, pi]
x_plot(:, 3) = mod(x_plot(:, 3) + pi, 2*pi) - pi;  % Wrap theta_g (ground heading)
x_plot(:, 6) = mod(x_plot(:, 6) + pi, 2*pi) - pi;  % Wrap theta_a (air heading)

x_nominal(:, 3) = mod(x_nominal(:, 3) + pi, 2*pi) - pi;  % Wrap theta_g (ground heading)
x_nominal(:, 6) = mod(x_nominal(:, 6) + pi, 2*pi) - pi;  % Wrap theta_a (air heading)

figure; % plot 1
for i = 1:size(x_plot, 2)
    subplot(size(x_plot, 2), 1, i);
    
    % Plot x_true in red
    plot(t, x_plot(:, i), 'r', 'LineWidth', 1.5); % 'r' specifies red color
    hold on;  % Keep the plot to overlay the next data set
    
    % Plot x_nominal in blue
    plot(t, x_nominal(:, i), 'b', 'LineWidth', 1.5); % 'b' specifies blue color
    
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
    
    % Add legend to clarify the plots
    if i == 1  % Add legend only in the first subplot to avoid repetition
        legend({'True State', 'Nominal State'}, 'Interpreter', 'latex');
    end
    
    hold off;  % Release the plot for next subplot operations
end
% Supertitle for all subplots
sgtitle('True State vs Nominal State Trajectories Comparison', 'Interpreter', 'latex');


% find total measurement vector
y_true = get_Y_noise(x_true,v_k');
y_nom = get_Y(x_nominal);

% Wrap the angles gamma_ag and gamma_ga to [-pi, pi]
y_true_plot = y_true;
y_nom_plot = y_nom;
y_true_plot(:, 1) = mod(y_true(:, 1) + pi, 2*pi) - pi;  % Wrap gamma_ag
y_true_plot(:, 3) = mod(y_true(:, 3) + pi, 2*pi) - pi;  % Wrap gamme_ga

y_nom_plot(:, 1) = wrapToPi(y_nom(:, 1));    % Wrap the first column (angles) of y_nom
y_nom_plot(:, 3) = wrapToPi(y_nom(:, 3));    % Wrap the third column (angles) of y_nom

% define ydata for testing
% ydata = ydata;    % default uses given ydata from canvas

figure;
for i = 1:size(y_true, 2)
    subplot(size(y_true, 2), 1, i);
    plot(t, y_true_plot(:, i), 'LineWidth', 1.5);
    hold on; % Keep the current plot
    plot(t, y_nom_plot(:, i), 'r', 'LineWidth', 1.5); % Red line for nominal data
    hold off; % Release the plot for new settings
    grid on;
    
    legend('True State', 'Nominal State', 'Location', 'best');

    % Adjust titles and labels based on the state variables
    if i == 1
        title('$\gamma_{ag}$ (Bearing UAV to UGV)', 'Interpreter', 'latex');
        ylabel('$\gamma_{ag}$ (rad)', 'Interpreter', 'latex');
    elseif i == 2
        title('$\rho_{ga}$ (Distance UGV to UAV)', 'Interpreter', 'latex');
        ylabel('$\rho_{ga}$ (m)', 'Interpreter', 'latex');
    elseif i == 3
        title('$\gamma_{ga}$ (Bearing UGV to UAV)', 'Interpreter', 'latex');
        ylabel('$\gamma_{ga}$ (rad)', 'Interpreter', 'latex');
    elseif i == 4
        title('$\xi_{a}$ (Easting of UAV)', 'Interpreter', 'latex');
        ylabel('$\xi_{a}$ (m)', 'Interpreter', 'latex');
    elseif i == 5
        title('$\eta_{a}$ (Northing of UAV)', 'Interpreter', 'latex');
        ylabel('$\eta_{a}$ (m)', 'Interpreter', 'latex');
    end
    xlabel('Time (s)', 'Interpreter', 'latex');
end
sgtitle('Truth Measurments vs. Time', 'Interpreter', 'latex');


%%%%%%%
% INITIALIZING EKF LOOP
%%%%%%%

%%%%%%%
% want to use a batch estimator to warm start the EKF



% okay... I don't think this function works
ydata = y_true';
dely = y_true - y_nom;

[del_x_0,P_plus_0] = get_init_conditions(dely,x_nominal,u_init,dt,Rtrue,10);
%nick's experiment:
x_init = x_init + del_x_0
P_init = P_plus_0

% my function
%delu = [0, 0, 0, 0]';
%[~, P_init] = warm_start(dely', x_nominal, delu, L, Rtrue);
%x_init = x_init + delx_0;

% human guess... doesn't result in good performance
%P_init = 0.1 * diag([1, 1, 1, 1, 1, 1]);

%%%%%%%
% gonna try and tune Q in here
% 0.6*Qtrue seems to perform the best
Qtrue = 0.6*Qtrue;
%{
Qtrue = [0.006, 0,      0,      0,      0,      0;
         0,     0.0006, 0,      0,      0,      0;
         0,     0,      0.0060, 0,      0,      0;
         0,     0,      0,      0.0006, 0,      0;
         0,     0,      0,      0,      0.0006, 0;
         0,     0,      0,      0,      0,      0.0060];
%%%%%%%
%}

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
    % Correct for wrapping issues in the first and third measurements (states)
    e_y(1, 1) = wrappedAngleDiff(ydata(1, k+1), y_min(1, 1));
    e_y(3, 1) = wrappedAngleDiff(ydata(3, k+1), y_min(3, 1));

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
