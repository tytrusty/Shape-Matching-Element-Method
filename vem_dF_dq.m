function dF_dq = vem_dF_dq(B, dM_dX)
    % Again, more work needed for quadratic deformation
    % dP_dq will then be in terms of X
    % dM_dX will be a vector of matrices
    % but I will dig this hole of laziness deeper
    n = size(B,1);      % # of points
    d = size(B,2); % dimension (2 or 3)
    
    dF_dq = zeros(d*d, d*n);
    
    for i = 1:n
        
        for j = 1:d
            dP_dq = zeros(d,n);
            dP_dq(j,:) = -(1/n);
            dP_dq(j,i) = (n-1)/n;
            dA_dq = dP_dq * B;
            dF_dq_i = dA_dq * dM_dX;
            dF_dq(:, d*(i-1) + j) = dF_dq_i(:);
        end
    end
end
