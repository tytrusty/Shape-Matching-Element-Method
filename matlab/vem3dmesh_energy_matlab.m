function [e, g, H] = vem3dmesh_energy_matlab(qdot_new, q, qdot, f_ext, x_fixed, ...
                                      vol, params, dF_dc, dF_dc_S, w_I, E, ...
                                      M, ME, L, P, J, ...
                                      k, n, d, x0_coms_size, k_stability, dt)
%%%%%                                    
% objective needs to return energy, gradient and hessian values
%%
   % Update position
    q_new = q + dt*qdot_new;
    x_new = reshape(P'*J*q_new,3,[]) + x_fixed;
    
    % Solve for polynomial coefficients (projection operators).
    c = vem3dmesh_polynomial_coefficients(x_new, L, E);
    
    neohookean_e =  vem3dmesh_neohookean_q(c, vol, params, dF_dc, dF_dc_S, d);
    
    PMP = P*M*P';
     
    e = 0.5*qdot_new'*J'*PMP*J*qdot_new - qdot_new'*J'*PMP*J*qdot + ...
        k_stability * x_new(:)' * ME * x_new(:) + ...
        0.5 * neohookean_e - ...
        qdot_new'*J'*f_ext; 

    if nargout > 1
        g_neohookean = vem3dmesh_neohookean_dq(x_new, c, vol, params, dF_dc, dF_dc_S, ME, L, ...
                                       k, n, d, x0_coms_size, k_stability);
                               
        g = J' * (PMP*J*(qdot_new - qdot) + ...
            + dt*P*g_neohookean + ...
            - f_ext);

        if nargout > 2  
            K = -vem3dmesh_neohookean_dq2(c, vol, params, dF_dc, w_I, k, n, ...
                                      x0_coms_size);
            K = L' * K * L;
          
            H = J' * (PMP + dt*dt*P*K*P') * J;
        end
    end
end

