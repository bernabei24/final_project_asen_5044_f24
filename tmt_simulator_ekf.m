% Extended Kalman Filter
clear;
clc;

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
[t, x_nominal] = ode45(dynamics_nominal, tvec, x_init); %perturbation not added


% Wrap the angles theta_g (x(3,:)) and theta_a (x(6,:)) to [-pi, pi]
%x_true(:, 3) = mod(x_true(:, 3) + pi, 2*pi) - pi;  % Wrap theta_g (ground heading)
%x_true(:, 6) = mod(x_true(:, 6) + pi, 2*pi) - pi;  % Wrap theta_a (air heading)

x_nominal(:, 3) = mod(x_nominal(:, 3) + pi, 2*pi) - pi;  % Wrap theta_g (ground heading)
x_nominal(:, 6) = mod(x_nominal(:, 6) + pi, 2*pi) - pi;  % Wrap theta_a (air heading)

% find total measurement vector
y_true = get_Y_noise(x_true,v_k');
y_nom = get_Y(x_nominal);
% define ydata for testing
% ydata = ydata;    % default uses given ydata from canvas
ydata = y_true';

%%%%%%%
% INITIALIZING EKF LOOP
%%%%%%%

P_init = 0.1 * diag([1, 1, 1, 1, 1, 1]);

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


%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%
%  NEES and NIS Testing Starts Here
%
%
%  TODO - Make this cleaner over the weekend.
%
%  
%  
%  This is a rough draft of the NEES and NIS testing
%  Just trying to get the logic down for now.
%
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

numberOfSimulationRuns = 20;
timeStepk = length(tvec)-1;  
NEESsshist = zeros(1, length(tvec));
NISsshist = zeros(1, length(tvec));


% Initialize the NEES and NIS samples
NEESsamps = zeros(numberOfSimulationRuns, length(tvec));
NISsamps = zeros(numberOfSimulationRuns, length(tvec));

% Loop through each simulation run
for ss = 1:numberOfSimulationRuns
    % Initialize the NEES and NIS history
    NEESsshist = zeros(1, length(tvec));
    NISsshist = zeros(1, length(tvec));

    %
    %  TODO - Figure out how to properly initialize TMT according to Dr. Nisar's instructions
    %
    %
    x_plus(:,1) = mvnrnd(x_plus(:,1), P_plus(:,:,1))';
    

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

    S = H * P_min(:, :, k+1) * H' + Rtrue;

    x_plus(:, k + 1) = x_min(:, k + 1) + K * e_y;
    P_plus(:, :, k + 1) = (eye(6) - K * H) * P_min(:, :, k+1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute the NEES statistic
    %
    % (NEES) error calculated here, that is,
    %
    %  e_x,k = (xhat - xtrue)'*P^+_k*(xhat - xtrue)
    %
    %
    NEESsshist(k) = (x_true(k,:)' - x_plus(:,k))'*P_plus(:,:,k+1)*(x_true(k,:)' - x_plus(:,k));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute the NIS statistic
    % 
    % (NIS) error calculated here, that is, 
    %
    %  e_y,k = (y - yhat)'*S^(-1)_k*(y - yhat)
    %
    NISsshist(k) = (e_y') / S * e_y;


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    
% Store the NEES and NIS samples for this simulation run
NEESsamps(ss, :) = NEESsshist;
NISsamps(ss, :) = NISsshist;

    


end

% Compute the mean and variance of the NEES and NIS samples
NEESmean = mean(NEESsamps, 1);
NISmean = mean(NISsamps, 1);

% Make alpha tuneable 
alphaNEES = 0.05;
alphaNIS = 0.05;

Nnx = numberOfSimulationRuns*length(F_k);
Nny = numberOfSimulationRuns*size(H,1);

%
%  N = number of simulation runs
%  n = degrees of freedom
%
% r1x and r2x are the lower and upper bounds of the confidence interval for the NEES statistic

% Compute the confidence intervals for the NEES and NIS statistics
r1x = chi2inv(alphaNEES/2, Nnx)./numberOfSimulationRuns;
r2x = chi2inv(1 - alphaNEES/2, Nnx)./numberOfSimulationRuns;

% r1y and r2y are the lower and upper bounds of the confidence interval for the NIS statistic
r1y = chi2inv(alphaNIS/2, Nny) ./ numberOfSimulationRuns;
r2y = chi2inv(1 - alphaNIS/2, Nny) ./ numberOfSimulationRuns;



figure(1)
plot(NEESmean,'ro','MarkerSize',6,'LineWidth',2),hold on
plot(r1x*ones(size(NEESmean)),'r--','LineWidth',2)
plot(r2x*ones(size(NEESmean)),'r--','LineWidth',2)
ylabel('NEES statistic, \bar{\epsilon}_x','FontSize',14)
xlabel('time step, k','FontSize',14)
xlim([1, length(tvec)-1])
title('NEES Estimation Results','FontSize',14)
legend('NEES @ time k', 'r_1 bound', 'r_2 bound')

figure(2)
plot(NISmean,'bo','MarkerSize',6,'LineWidth',2),hold on
plot(r1y*ones(size(NISmean)),'b--','LineWidth',2)
plot(r2y*ones(size(NISmean)),'b--','LineWidth',2)
ylabel('NIS statistic, \bar{\epsilon}_y','FontSize',14)
xlabel('time step, k','FontSize',14)
xlim([1, length(tvec)-1])
title('NIS Estimation Results','FontSize',14)
legend('NIS @ time k', 'r_1 bound', 'r_2 bound')
