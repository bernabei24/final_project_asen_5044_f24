function [Y] = get_Y_noise(x, noise)
% Inputs
    % x: 1001x6 matrix of states
    % noise: 1001x5 matrix of measurement noise
% Outputs
    % Y: 1001x5 matrix of sensor data with noise added
    Y = [];
    for k = 1:size(x, 1) 
        x_k = x(k, :)';
        y_k = [atan2((x_k(5) - x_k(2)), (x_k(4) - x_k(1))) - x_k(3);
               sqrt((x_k(1) - x_k(4))^2 + (x_k(2) - x_k(5))^2);
               atan2((x_k(2) - x_k(5)), (x_k(1) - x_k(4))) - x_k(6);
               x_k(4);
               x_k(5)];
        % Add noise to the output y_k
        y_k = y_k + noise(k, :)';
        Y = [Y; y_k'];
    end
end
