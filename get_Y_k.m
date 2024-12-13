function Y_k = get_Y_k(x_k)
% inputs
    % x_k: 6-element state vector at time k
% outputs
    % Y_k: 5-element sensr values at time k
    Y_k = [atan2((x_k(5) - x_k(2)), (x_k(4) - x_k(1))) - x_k(3);
           sqrt((x_k(1) - x_k(4))^2 + (x_k(2) - x_k(5))^2);
           atan2((x_k(2) - x_k(5)), (x_k(1) - x_k(4))) - x_k(6);
           x_k(4);
           x_k(5)];
end

