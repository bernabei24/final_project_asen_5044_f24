function [B_bar] = get_Bbar(u,x,L)
%inputs: 
    %u: vector of inputs in form u = [ug ua]^T = [vg phig va phia]^T
    %x: vector of inputs in form x = [Eg Ng thetag Ea Na thetaa]^T
    %L: length between front wheels of ground vehicle
    B_bar = [cos(x(3,1)) 0 0 0;
             sin(x(3,1)) 0 0 0;
             tan(u(2,1))/L (u(1,1)/L)*( sec(u(2,1))^2 ) 0 0;
             0 0 cos(x(6,1)) 0;
             0 0 sin(x(6,1)) 0;
             0 0 0 1];
    
end