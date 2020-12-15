function J = vem_jacobian(B, Q, n, d, N, E)
    J = zeros(d,N*d, size(Q,2));
    
    dcom_dq = repelem(-1/N,n+1);
    for i = 1:size(Q,2)
        BM_i = [B * Q(:,i); -1]; % Append 1 to BM
        dcom_i = dcom_dq * BM_i;
        for j = 1:d
           J(j,j:d:end,i)  = dcom_i;
        end
        
        for j = 1:numel(E)
            for k = 1:d
                idx = d*(E(j)-1);
                J(k, idx + k, i) = squeeze(J(k, idx + k, i)) + BM_i(j);
            end
        end
    end
end

