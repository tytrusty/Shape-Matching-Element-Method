function dM_dX = monomial_basis_grad(x, x_com, order)
    d = size(x,1); % dimension
    m = size(x,2);
    
    % Compute number of monomials
    k = basis_size(d, order);
    
    dM_dX = zeros(m, k, d);
    Q = x - x_com;
    
    % TODO -- make this general
    for ii = 1:m
        dMi_dX = zeros(k,d);
        dMi_dX(1:d,:) = eye(d);
        if order == 2
            if d==2
                dMi_dX(3,:) = [2*Q(1,ii) 0];
                dMi_dX(4,:) = [0 2*Q(2,ii)];
                dMi_dX(5,:) = [Q(2,ii) Q(1,ii)];
            else
                dMi_dX(4,:) = [2*Q(1,ii) 0 0];
                dMi_dX(5,:) = [0 2*Q(2,ii) 0];
                dMi_dX(6,:) = [0 0 2*Q(3,ii)];
                dMi_dX(7,:) = [Q(2,ii) Q(1,ii) 0];
                dMi_dX(8,:) = [0 Q(3,ii) Q(2,ii)];
                dMi_dX(9,:) = [Q(3,ii) 0 Q(1,ii)];
            end
        end
        dM_dX(ii,:,:) = dMi_dX;
    end
end

