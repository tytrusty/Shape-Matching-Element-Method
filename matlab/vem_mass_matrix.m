function M = vem_mass_matrix2(Y, W, W_S, L)
    m = size(Y,1);
    M = zeros(size(L,1), size(L,1));
    d = size(Y,2);
    for i=1:m
    	Yi = squeeze(Y(i,:,:))*W{i} * W_S{i}; % weighed monomial basis
        Yi(:, end-d+1:end) = eye(d);
        M = M +  Yi'*Yi;
    end
    M = L' * M * L;

end