function [f] = vem3dmesh_neohookean_dq_matlab(x, c, vol, params, dF_dc, dF_dc_S, ME, L, ...
                                       k, n, d, m, x0_coms_size, k_stability)
                                    
   % Force vector
    dV_dq = zeros(d*(k*n + x0_coms_size),1);

    % Computing force dV/dq for each point.
    % TODO: move this to C++ :)
    for i = 1:m
        % Deformation Gradient
        F = dF_dc{i} * dF_dc_S{i} * c;
        F = reshape(F,d,d);

        % Force vector
        dV_dF = neohookean_tet_dF(F, params(i,1), params(i,2));
        dV_dq = dV_dq +  dF_dc_S{i}' * dF_dc{i}' * dV_dF * vol(i);
    end
    dV_dq = L' * dV_dq;

    % Error correction force
    f_error = - 2 * ME * x(:);
    f_error = k_stability*f_error(:);

    % Force from potential energy.
    f_internal = -dV_dq;

    f = f_internal + f_error;
end

