function M = vem_error_matrix(B, Q, n, d)
    M = zeros(n*d,n*d);
    dP_dq = zeros(d*n,d*(n+1));
    
    for i = 1:n
        for j = 1:d
            idx1 = (i-1)*d + j;
            idx2 = (j-1)*(n+1) + 1; 
            dP_dq(idx1,idx2:idx2+(n-1)) = 1/n;
            dP_dq(idx1,idx2+(i-1)) = -(n-1)/n;
            
            % center of mass contribution
            dP_dq(idx1,idx2+n) = -1/n;
        end
    end
    for i = 1:size(Q,2)
        BM = [B * Q(:,i); 1]; % Append 1 to BM
        zero = zeros(size(BM));
        BM = [BM zero; zero BM];
        
        % Build new dP_dq mat
        dP = dP_dq;
        idx1 = (i-1)*d;
        dP(idx1 + 1, n+1)= (n+1)/n;
        dP(idx1 + 2, 2*(n+1))= (n+1)/n;
        
        Mi = dP * BM;
        %Mi = dP_dq * BM;
        M = M + Mi*Mi';
    end
end

