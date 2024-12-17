function [x_0, P_0] = warm_start(y, x_nom, R)

    y = y(:, 1:10);
    y = reshape(y, [50, 1]);
    R = blkdiag(R, R, R, R, R, R, R, R, R, R);

    H = [];
    for k = 1:10
        C_nom_k = get_Cbar(x_nom(k, :)');
        H_k = C_nom_k;

        H = [H; H_k];
    end

    P_0 = (H' * R^(-1) * H)^(-1)
    x_0 = P_0 * H' * R^(-1) * y

end

