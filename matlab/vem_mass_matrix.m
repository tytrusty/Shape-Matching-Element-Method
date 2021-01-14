function M = vem_mass_matrix(Y, W, W_S, L, mass)
    m = size(Y,1);
    M = zeros(size(L,1), size(L,1));
    d = size(Y,2);
    for i=1:m
%     	Yi = squeeze(Y(i,:,:))*W{i} * W_S{i}; % weighed monomial basis
%         Yi(:, end-d+1:end) = eye(d);
%         M = M +  Yi'*Yi * mass(i);
        M = M +  Y{i}'*Y{i} * mass(i);
    end
    M = L' * M * L;
end