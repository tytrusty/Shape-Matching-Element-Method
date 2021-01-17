function ME = vem_error_matrix(Y, Y_S, L, d)
    m = size(Y,1);    
    M1 = zeros(size(L,1),size(L,1));
    M2 = zeros(size(L,1),size(L,2));
    M3 = eye(size(L,2));
    tic
    for i=1:m
        M1 = M1 + Y_S{i}'*(Y{i}'*Y{i})*Y_S{i};
        M2(:,d*(i-1)+1:d*i) = Y_S{i}'*Y{i}';
    end
    ME = L'*M1*L - L'*M2 - M2'*L + M3;
end