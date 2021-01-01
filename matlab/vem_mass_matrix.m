function M = vem_mass_matrix(B, Q, w, d, N, E)
    m = size(Q,2);
    M = zeros(d*(N+1), d*(N+1));
    J = zeros(d,d*(N+1),m);
    for i = 1:size(E,1)
        w_i = reshape(w(:,i), [1 1 m]);
        J = J + bsxfun(@times, vem_jacobian(B{i},Q,d,N,E{i}), w_i);
    end
    
    for j=1:m        
        M = M + J(:,:,j)'*J(:,:,j);

    end
end