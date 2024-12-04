function [A_bar] = get_Abar(u,x)
%inputs: 
    %u: vector of inputs in form u = [ug ua]^T = [vg phig va phia]^T
    %x: vector of inputs in form x = [Eg Ng thetag Ea Na thetaa]^T

    A_bar = [0 0 -u(1,1)*sin(x(3,1)) 0 0 0;
             0 0  u(1,1)*cos(x(3,1)) 0 0 0;
             0 0  0                  0 0 0;
             0 0 0 0 0  u(3,1)*sin(x(6,1));
             0 0 0 0 0 -u(3,1)*cos(x(6,1));
             0 0 0 0 0 0];

end