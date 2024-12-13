%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% .__   __.  _______  _______     _______.   .__   __.  __       _______.        _______.  ______ .______       __  .______   .___________.
% |  \ |  | |   ____||   ____|   /       |   |  \ |  | |  |     /       |       /       | /      ||   _  \     |  | |   _  \  |           |
% |   \|  | |  |__   |  |__     |   (----`   |   \|  | |  |    |   (----`      |   (----`|  ,----'|  |_)  |    |  | |  |_)  | `---|  |----`
% |  . `  | |   __|  |   __|     \   \       |  . `  | |  |     \   \           \   \    |  |     |      /     |  | |   ___/      |  |     
% |  |\   | |  |____ |  |____.----)   |      |  |\   | |  | .----)   |      .----)   |   |  `----.|  |\  \----.|  | |  |          |  |     
% |__| \__| |_______||_______|_______/       |__| \__| |__| |_______/       |_______/     \______|| _| `._____||__| | _|          |__|     
% 
% 
% This script is attempting to implement the suggested model for NEES and NIS testing that Dr. Ahmed gave us in class.
% Most of these components have been implemented by Whit and Nick, just trying to ensure all entities in 
% Dr. Ahmed's suggested model are present.
%
%  See lecture 29 slide 14 for his suggested model
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% .___________..___  ___. .___________.        _______. __  .___  ___.  __    __   __          ___   .___________.  ______   .______      
% |           ||   \/   | |           |       /       ||  | |   \/   | |  |  |  | |  |        /   \  |           | /  __  \  |   _  \     
% `---|  |----`|  \  /  | `---|  |----`      |   (----`|  | |  \  /  | |  |  |  | |  |       /  ^  \ `---|  |----`|  |  |  | |  |_)  |    
%     |  |     |  |\/|  |     |  |            \   \    |  | |  |\/|  | |  |  |  | |  |      /  /_\  \    |  |     |  |  |  | |      /     
%     |  |     |  |  |  |     |  |        .----)   |   |  | |  |  |  | |  `--'  | |  `----./  _____  \   |  |     |  `--'  | |  |\  \----.
%     |__|     |__|  |__|     |__|        |_______/    |__| |__|  |__|  \______/  |_______/__/     \__\  |__|      \______/  | _| `._____|
%                                                                                                                                        
% 
% The Truth Model Tester (TMT) has been handled by the work done with ODE45 simulation in the previous scripts, I believe?  These simulation scripts 
% provide the needed y^k coming from the Truth Model (TM) Simulator.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  __  ___  _______      ______   ______    _______   _______ 
% |  |/  / |   ____|    /      | /  __  \  |       \ |   ____|
% |  '  /  |  |__      |  ,----'|  |  |  | |  .--.  ||  |__   
% |    <   |   __|     |  |     |  |  |  | |  |  |  ||   __|  
% |  .  \  |  |        |  `----.|  `--'  | |  '--'  ||  |____ 
% |__|\__\ |__|         \______| \______/  |_______/ |_______|
%
% Again this is implmented already in the previous scripts, I believe?  The Part4Wrapper.m houses this logic
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% .__   __.  _______  _______     _______.     ___.__   __.  __       _______.        ___       __        _______   ______   
% |  \ |  | |   ____||   ____|   /       |    /  /|  \ |  | |  |     /       |       /   \     |  |      /  _____| /  __  \  
% |   \|  | |  |__   |  |__     |   (----`   /  / |   \|  | |  |    |   (----`      /  ^  \    |  |     |  |  __  |  |  |  | 
% |  . `  | |   __|  |   __|     \   \      /  /  |  . `  | |  |     \   \         /  /_\  \   |  |     |  | |_ | |  |  |  | 
% |  |\   | |  |____ |  |____.----)   |    /  /   |  |\   | |  | .----)   |       /  _____  \  |  `----.|  |__| | |  `--'  | 
% |__| \__| |_______||_______|_______/    /__/    |__| \__| |__| |_______/       /__/     \__\ |_______| \______|  \______/  
%
%  
%  Algorithm for NEES and NIS testing:
%  1.  Start Part4Wrapper.m and ensure that we are getting y^k from the Truth Model (TM) Simulator and ground truth x^k
%  2.  For k = 1,2,...,N
%  3.  Compute the Kalman gain
%  4.  Update the state estimate
%  5.  Update the state covariance
%  6.  Compute the NEES and NIS statistics (we are computing epsilon_x,k and epsilon_y,k in each time step)
%  7.  Repeat
%
%  The NEES and NIS statistics are computed as follows:
%  NEES = (xhat - xtrue)'*P^(-1)*(xhat - xtrue)
%  NIS = (y - yhat)'*R^(-1)*(y - yhat)
%
%  Where xhat is the state estimate, xtrue is the true state, P is the state covariance, y is the measurement, yhat is the predicted measurement, and S is the measurement covariance.
%
%  The NEES and NIS statistics are then used to determine if the filter is consistent.  If the filter is consistent, the NEES and NIS statistics should be distributed according to a chi-squared distribution with the appropriate degrees of freedom.
%


%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


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
x_init = [Eg_init Ng_init thetag_init Ea_init Na_init thetaa_init]';
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
x_true(:, 3) = mod(x_true(:, 3) + pi, 2*pi) - pi;  % Wrap theta_g (ground heading)
x_true(:, 6) = mod(x_true(:, 6) + pi, 2*pi) - pi;  % Wrap theta_a (air heading)

x_nominal(:, 3) = mod(x_nominal(:, 3) + pi, 2*pi) - pi;  % Wrap theta_g (ground heading)
x_nominal(:, 6) = mod(x_nominal(:, 6) + pi, 2*pi) - pi;  % Wrap theta_a (air heading)


% find total measurement vector
y_true = get_Y_noise(x_true,v_k');
y_nom = get_Y(x_nominal);

% Wrap the angles gamma_ag and gamma_ga to [-pi, pi]
y_true(:, 1) = mod(y_true(:, 1) + pi, 2*pi) - pi;  % Wrap gamma_ag
y_true(:, 3) = mod(y_true(:, 3) + pi, 2*pi) - pi;  % Wrap gamme_ga


%INITIALIZING FOR LKF LOOP
%pretty sure this is what gamme looks like although not positive
gamma = eye(size(x_true,2));
omega_matrix = dt * gamma;

dely = y_true - y_nom;
delx_plus = zeros(6,1001);
delx_minus = zeros(6,1001);
P_plus = zeros(6,6,1001);
P_minus = zeros(6,6,1001);
K = zeros(6,5,1001);

delx_plus(:,1) = zeros(6,1); %adjustable
P_plus(:,:,1) = zeros(6,6); %adjustable
Q_filter = Qtrue;


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
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %  This is done by the Part4Wrapper.m script (LKF Implementation), 
    %  the Kalman call will be made into a function call.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Perform the prediction step
    %Compute jacobians at each time step from nominal trajectory
    for k = 1:timeStepk %k=1 represents t = 0

        % Use CT jacobians to find A and B at t = t_k and nom[k]
        A_nom_k = get_Abar(u_init, x_nominal(k, :)');
        B_nom_k = get_Bbar(u_init, x_nominal(k, :)', L);
        C_nom_k = get_Cbar(x_nominal(k+1, :)'); %need H_k+1 not H_k
        
        % Eulerized estimate of DT Jacobians
        F_k = dt * A_nom_k;
        G_k = dt * B_nom_k;
        H_k_plus_1 = C_nom_k;%need H_k+1 not H_k
        
        %LINEAR KALMAN FILTER

        %time update/prediction step for time k+1
        delx_minus(:,k+1) = F_k * delx_plus(:,k) + G_k * [0,0,0,0]';
        P_minus(:,:,k+1) = F_k*P_plus(:,:,k)*F_k' + omega_matrix*Qtrue*omega_matrix';

        innov_kp1 = dely(k+1,:)'; 

        S = H_k_plus_1 * P_minus(:,:,k+1) * H_k_plus_1' + Rtrue;

        %measurement update/correction step for time k+1
        K(:,:,k+1) = P_minus(:,:,k+1) * H_k_plus_1' * (S)^(-1);
        
        delx_plus(:,k+1) = delx_minus(:,k+1) + K(:,:,k+1) * ...
                            (dely(k+1,:)' - H_k_plus_1 * delx_minus(:,k+1) );
        
        P_plus(:,:,k+1) = (eye(6) - K(:,:,k+1) * H_k_plus_1) * P_minus(:,:,k+1);


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Compute the NEES statistic
        %
        % (NEES) error calculated here, that is,
        %
        %  e_x,k = (xhat - xtrue)'*P^+_k*(xhat - xtrue)
        %
        %
        NEESsshist(k) = (x_true(k,:)' - delx_plus(:,k))'*P_plus(:,:,k+1)*(x_true(k,:)' - delx_plus(:,k));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Compute the NIS statistic
        % 
        % (NIS) error calculated here, that is, 
        %
        %  e_y,k = (y - yhat)'*S^(-1)_k*(y - yhat)
        %
        NISsshist(k) = (innov_kp1') / S * innov_kp1;


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
Nny = numberOfSimulationRuns*size(H_k_plus_1,1);

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
















