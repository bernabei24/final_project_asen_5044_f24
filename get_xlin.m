function [xlin] = get_xlin(F,G,x_init,u_init,timesteps)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
    xlin = zeros(size(x_init,1),timesteps);
    for k = 1:timesteps
        if k == 1
            xlin(:,1) = F*x_init + G*u_init;
        else
            xlin(:,k) = F*xlin(:,k-1) + G*u_init;
        end
    end


end