function M = vem_mass_matrix(Y, Y_S, L, mass)
    m = size(Y,1);
    M = zeros(size(L,1), size(L,1));
    for i=1:m
        M = M + Y_S{i}'*(Y{i}'*Y{i})*Y_S{i} * mass(i);
    end
    M = L' * M * L;
end