function [Y] = get_Y(x)
% inputs
    % x: 1001x6 matrix of states
% outputs
    % Y: 1001x5 matrix of sensor data
    Y = [];
    for k = 1:length(x)
        x_k = x(k, :)';
        y_k = [atan2((x_k(5) - x_k(2)), (x_k(4) - x_k(1))) - x_k(3);
               sqrt((x_k(1) - x_k(4))^2 + (x_k(2) - x_k(5))^2);
               atan2((x_k(2) - x_k(5)), (x_k(1) - x_k(4))) - x_k(6);
               x_k(4);
               x_k(5)];
        Y = [Y;
             y_k'];
    end
end

