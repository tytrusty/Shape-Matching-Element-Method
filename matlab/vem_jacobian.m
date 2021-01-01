function J = vem_jacobian(B, Q, d, N, E)
    J = zeros(d,d*(N+1), size(Q,2));
    
    for i = 1:size(Q,2)
        BM_i = B * Q(:,i);

        for j = 1:numel(E)
            for k = 1:d
                idx = d*(E(j)-1);
                J(k, idx + k, i) = squeeze(J(k, idx + k, i)) + BM_i(j);
                J(k, d*N + k, i) = squeeze(J(k, d*N + k, i)) - BM_i(j);
            end
        end
        J(:,d*N+1:d*(N+1),i) = J(:,d*N+1:d*(N+1),i) + eye(d);
        
    end
end

