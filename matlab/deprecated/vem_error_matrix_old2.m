function ME = vem_error_matrix(B, Q, w, d, N, E)
    ME = zeros(d*N, d*N);
    for i = 1:size(E,1)
        n=size(B{i},1);

        % Stability term
        J = vem_jacobian(B{i},Q,n,d,N,E{i});
        for j=1:N
            I = zeros(d,d*N);
            I(:, d*j-(d-1):d*j) = eye(d);
            
            if w(j,i) > 1e-12
                Jtmp = w(j,i) * (I - J(:,:,j));
                ME = ME + Jtmp'*Jtmp;
            end
        end
        
    end
end