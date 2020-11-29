function J = vem_jacobian(B, Q, n, d, nn)
    if nargin <5
        nn=n
    end
    J = zeros(d,n*d, size(Q,2));
    
    %M = zeros(n*d, n*d);
    
    dP_dq = zeros(d,d*(n+1)*n);
    
    % See VEM_A_func
    for i = 1:n
        for j = 1:d
            idx = (i-1)*d*(n+1) + (j-1)*(n+1) + 1;
            dP_dq(j, idx:idx+(n-1)) = -1/nn;
            dP_dq(j, idx+(i-1)) = (nn-1)/nn;
            dP_dq(j, idx+n) = 1/nn;
        end
    end
    
    for i = 1:size(Q,2)
        BM = zeros(d*n*(n+1),d*n); 
        BM_i = [B * Q(:,i); 1]; % Append 1 to BM
        for j = 1:n*d
            idx = (j-1)*(n+1) + 1;
           BM(idx:idx+n,j)=BM_i;
        end
        
        J(:,:,i) = dP_dq * BM;
        %M = M + J'*J;
    end
end

