function ME = vem_error_matrix2(Y, W, W_S, L)
    m = size(Y,1);
    ME = zeros(size(L,2), size(L,2));
    d = size(Y,2);

    for i=1:m
    	Yi = squeeze(Y(i,:,:))*W{i} * W_S{i}; % weighed monomial basis
        Yi(:, end-d+1:end) = eye(d);
        
        I = zeros(d,size(L,2));
        I(:, d*(i-1)+1: d*i) = eye(d);
        
        J = I - Yi*L;
        
        ME = ME + J'*J;
    end     
end