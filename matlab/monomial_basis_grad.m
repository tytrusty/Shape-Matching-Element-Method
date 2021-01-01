function dM_dX = monomial_basis_grad(x, x_com, order)
    d = size(x,1); % dimension
    m = size(x,2);
    
    % Compute number of monomials
    k=0;
    for i=1:order
        k=k+nchoosek(d+i-1,i);
    end
    
    dM_dX = zeros(m,k,d);
    Q = x - x_com;
    
    % TODO -- make this general
    for i = 1:m
        dMi_dX(1:d,:) = eye(d);
        if order == 2
            if d==2
                dMi_dX(3,:) = [2*Q(1,i) 0];
                dMi_dX(4,:) = [0 2*Q(2,i)];
                dMi_dX(5,:) = [Q(2,i) Q(1,i)];
            else
                dMi_dX(4,:) = [2*Q(1,i) 0 0];
                dMi_dX(5,:) = [0 2*Q(2,i) 0];
                dMi_dX(6,:) = [0 0 2*Q(3,i)];
                dMi_dX(7,:) = [Q(2,i) Q(1,i) 0];
                dMi_dX(8,:) = [0 Q(3,i) Q(2,i)];
                dMi_dX(9,:) = [Q(3,i) 0 Q(1,i)];
            end
        end
        dM_dX(i,:,:) = dMi_dX;
    end
end

