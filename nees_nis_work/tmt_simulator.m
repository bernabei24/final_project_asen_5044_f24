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

% Function to compute the NEES and NIS statistics
function [NEESsamps, NISsamps] = NEESNIS(NEESsshist, NISsshist, numberOfSimulationRuns, timeStepk, xk_truehist, mk_filt_hist, Pk_filt_hist, ykhist, H, Rtrue, Qkf)
    
    % Initialize the NEES and NIS samples
    NEESsamps = zeros(numberOfSimulationRuns, length(timeStepk));
    NISsamps = zeros(numberOfSimulationRuns, length(timeStepk));
    
    % Loop through each simulation run
    for ss = 1:numberOfSimulationRuns
        % Initialize the NEES and NIS history
        NEESsshist = zeros(1, length(timeStepk));
        NISsshist = zeros(1, length(timeStepk));
        
        % Loop through each time step
        for k = 1:length(timeStepk)


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            %  This is done by the Part4Wrapper.m script (LKF Implementation), 
            %  trying to figure out if we want to integrate in the Wrapper or make this a function call
            %
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Perform the prediction step
            mkp1_minus = timeStepk*mk_filt_hist(:,k) + G*u(jj);
            Pkp1_minus = timeStepk*Pk_filt_hist(:,:,k)*timeStepk' + Qkf;
            
            % Compute the Kalman gain
            Pyykp1 = H*Pkp1_minus*H' + Rtrue;
            Pyykp1 = 0.5*(Pyykp1 + Pyykp1');
            Kkp1 = Pkp1_minus*H'/(H*Pkp1_minus*H' + Rtrue);
            
            % Perform the measurement update step
            ykp1_report = ykhist(:,jj);
            ykp1_pred = H*mkp1_minus;
            innov_kp1 = ykp1_report - ykp1_pred;
            mkp1_plus = mkp1_minus + Kkp1*innov_kp1;
            Pkp1_plus = (eye(2) - Kkp1*H)*Pkp1_minus;
            
            % Update the state estimate and covariance
            mk_filt_hist(:,k) = mkp1_plus;
            Pk_filt_hist(:,:,k) = Pkp1_plus;
            
 
            invPkp1 = inv(Pkp1_plus);
            invPyykp1 = inv(Pyykp1);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            
            %
            % (NEES) error calculated here, that is,
            %
            %  e_x,k = (xhat - xtrue)'*P^+_k*(xhat - xtrue)
            %
            %
            NEESsshist(k) = (xk_truehist(:,k) - mkp1_plus)'*invPkp1*(xk_truehist(:,k) - mkp1_plus);
            
            % 
            % (NIS) error calculated here, that is, 
            %
            %  e_y,k = (y - yhat)'*R^(-1)_k*(y - yhat)
            %
            NISsshist(k) = innov_kp1'*invPyykp1*innov_kp1;


        end

        % Store the NEES and NIS samples
        NEESsamps(ss, :) = NEESsshist;
        NISsamps(ss, :) = NISsshist;


    end

end

% DEAD CODE, JUST WORKING OUT DESIGN AT THIS POINT, DON'T RUN THIS CODE
% Get NEESsamp and NISsamp
[NEESsamps, NISsamps] = NEESNIS(NEESsshist, NISsshist, numberOfSimulationRuns, timeStepk, xk_truehist, mk_filt_hist, Pk_filt_hist, ykhist, H, Rtrue, Qkf);

% Compute the mean and variance of the NEES and NIS samples
NEESmean = mean(NEESsamps, 1);
NISmean = mean(NISsamps, 1);

% Make alpha tuneable 
alphaNEES = 0.05;
alphaNIS = 0.05;


% Do the nulkl hypothesis testing


%
%  N = number of simulation runs
%  n = total number of time step's k
%
% r1x and r2x are the lower and upper bounds of the confidence interval for the NEES statistic

% Compute the confidence intervals for the NEES and NIS statistics
r1x = chi2inv(alphaNEES/2, numberOfSimulationRuns*timeStepk) ./ numberOfSimulationRuns;
r2x = chi2inv(1 - alphaNEES/2, numberOfSimulationRuns*timeStepk) ./ numberOfSimulationRuns;

% r1y and r2y are the lower and upper bounds of the confidence interval for the NIS statistic
r1y = chi2inv(alphaNIS/2, numberOfSimulationRuns*timeStepk) ./ numberOfSimulationRuns;
r2y = chi2inv(1 - alphaNIS/2, numberOfSimulationRuns*timeStepk) ./ numberOfSimulationRuns;

% TODO - Plot the NEES and NIS statistics with the confidence intervals











