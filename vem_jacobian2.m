function J = vem_jacobian2(B, Q, n, d, N, E)
    J = zeros(d,N*d, size(Q,2));
    
    %M = zeros(n*d, n*d);
    
    dP_dq = zeros(d,d*(n+1)*N);
    
    % See VEM_A_func
    for i = 1:N
        for j = 1:d
            idx = (i-1)*d*(n+1) + (j-1)*(n+1) + 1;
            dP_dq(j, idx:idx+(n-1)) = -1/N;
%             dP_dq(j, idx+(i-1)) = (N-1)/N;
            dP_dq(j, idx+n) = 1/N;
        end
    end
    
    for i = 1:numel(E)
        for j = 1:d
            idx = (E(i)-1)*d*(n+1) + (j-1)*(n+1) + 1;
            dP_dq(j, idx+(i-1)) = (N-1)/N;
        end
        
    end
    
    for i = 1:size(Q,2)
        BM = zeros(d*N*(n+1),d*N); 
        BM_i = [B * Q(:,i); 1]; % Append 1 to BM
        for j = 1:N*d
            idx = (j-1)*(n+1) + 1;
           BM(idx:idx+n,j)=BM_i;
        end
        
        J(:,:,i) = dP_dq * BM;
        %M = M + J'*J;
    end
end

