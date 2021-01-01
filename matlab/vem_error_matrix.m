function ME = vem_error_matrix(B, Q, w, d, N, E)
    ME = zeros(d*(N+1), d*(N+1));
    
    J = zeros(d,d*(N+1),N); % should this final dimension be N?
    for i = 1:size(E,1)
        w_i = reshape(w(:,i), [1 1 N]);
        J = J + bsxfun(@times, vem_jacobian(B{i},Q,d,N,E{i}), w_i);
    end

    for j=1:N
        I = zeros(d,d*(N+1));
        I(:, d*j-(d-1):d*j) = eye(d);
        Jtmp = I - J(:,:,j);
        ME = ME + Jtmp'*Jtmp;
    end
        
end