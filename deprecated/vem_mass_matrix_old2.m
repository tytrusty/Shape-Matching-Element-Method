function M = vem_mass_matrix(B, Q, w, d, N, E)
    M = zeros(d*N, d*N);
    for i = 1:size(E,1)
        n=size(B{i},1);

        J = vem_jacobian(B{i},Q,n,d,N,E{i});
        for j=1:size(Q,2)
            if w(j,i) > 1e-4
                Jtmp = w(j,i) * J(:,:,j);
                M = M + Jtmp'*Jtmp;
            end
        end
    end
end