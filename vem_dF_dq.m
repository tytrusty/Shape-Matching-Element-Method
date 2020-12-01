function dF_dq = vem_dF_dq(B, dM_dX)
    n = size(B,1);      % # of points
    d = size(dM_dX,2);  % dimension (2 or 3)
    
    dF_dq = zeros(d*d, d*n);
    for i = 1:n
        
        for j = 1:d
            dP_dq = zeros(d,n);
            dP_dq(j,:) = -(1/n);
            dP_dq(j,i) = (n-1)/n;
%             dP_dq(j,i) = 1;
            dA_dq = dP_dq * B;
            dF_dq_i = dA_dq * dM_dX;
            dF_dq(:, d*(i-1) + j) = dF_dq_i(:);
        end
    end
end
