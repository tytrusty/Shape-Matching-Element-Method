function M = vem_mass_matrix(B, Q, w, d, N, E)
    m = size(Q,2);
    M = zeros(d*N, d*N);
    J = zeros(d,d*N,m);
    for i = 1:size(E,1)
        n=size(B{i},1);
        w_i = reshape(w(:,i), [1 1 m]);
        J = J + bsxfun(@times, vem_jacobian(B{i},Q,n,d,N,E{i}), w_i);
    end
    for j=1:m
    	M = M + J(:,:,j)'*J(:,:,j);
    end
end