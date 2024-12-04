function [F,G] = get_FG(A_bar,B_bar, dt)
%inputs:    
    %A_bar matrix
    %B_bar matrix
    %dt: time interval in seconds
%outputs
    Ahat = [A_bar B_bar; zeros(4,10)];
    ZZ = expm(Ahat*dt);
    F = ZZ(1:6,1:6);
    G = ZZ(1:6,7:10);
    
end