function dP_dq = vem_dP_dq(n,d)
    %n = size(B,1);      % # of points
    %d = size(dM_dX,2);  % dimension (2 or 3)
    
    dP_dq = zeros(d*d*n,n);
    for i = 1:n
        for j = 1:d
            idx = d*d*i + d*j - d*d - 1 + (j-1);
            dP_dq(idx,:) = -(1/n);
            dP_dq(idx,i) = (n-1)/n;
        end
    end
end
