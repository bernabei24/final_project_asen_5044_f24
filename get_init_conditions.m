function [del_x_0,P_plus_0] = get_init_conditions(dely,x_nominal,u_init,dt,Rtrue,meas_considered)
    %dely — measurements to cpmute delx_0
    %x_nominal — needed to compute F matrixes and H matrixes
    %u_init — needed to compute F matrixes and H matrixes
    %dt — needed to compute F matrixes and H matrixes
    %Rtrue — used to create Raug matrix
    %meas_considered — number of dely measurements used to compute x_0
    
    p = size(dely,2); %number of variables in measurement vector
    n = size(x_nominal,2); %number of variables in state vector
    
    
    big_dely = zeros(meas_considered*p,1);
    big_alpha = zeros(meas_considered*p,n);
    big_R = zeros(meas_considered*p,meas_considered*p);
    STM = eye(6); %help build big_alpha
    
    for k = 1:meas_considered
        
        A_nom_k = get_Abar(u_init, x_nominal(k, :)');
        C_nom_k = get_Cbar(x_nominal(k, :)');

        F_k = eye(6) + dt * A_nom_k;
        H_k = C_nom_k;
        
        big_dely((5*k-4):(5*k),:) = dely(k,:)';
        if k == 1
            big_alpha((5*k-4):(5*k),:) = H_k;
            big_R( (5*k-4):(5*k),(5*k-4):(5*k) ) = Rtrue;
            STM = F_k * STM;
        else
            big_alpha((5*k-4):(5*k),:) = H_k * STM;
            STM = F_k * STM;
            big_R( (5*k-4):(5*k),(5*k-4):(5*k) ) = Rtrue;
        end
    end

    del_x_0 = ( ( big_alpha' * big_R^-1 * big_alpha)^-1 *big_alpha' * big_R^-1 ) * big_dely;
    P_plus_0 = (big_alpha' * big_R^-1 * big_alpha)^-1;

            
end