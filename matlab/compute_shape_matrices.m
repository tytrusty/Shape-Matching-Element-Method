function [B,L,LM] = compute_shape_matrices(x0, x0_com, E, order)
    d = size(x0,1); % dimension

    % Compute number of monomials
    k = basis_size(d, order);

    Q = cell(numel(E),1);
    B = cell(numel(E),1);

    M_rows = 0;
    for i=1:numel(E)
        x = x0(:,E{i});
        Q{i} = monomial_basis(x, x0_com, order);
        M_rows = M_rows + numel(E{i});
    end
    
    M_cols = d*k*numel(E) + d;
    M = zeros(d * M_rows, M_cols);

    col_idx = 0;
    row_idx = 0;
    for i=1:numel(E)
        for j=1:numel(E{i})
            for n=1:d
                col_start = col_idx + (n-1)*k + 1;
                col_end   = col_idx + n*k;
                M(row_idx+n, col_start:col_end) = Q{i}(:,j);
            end
            M(row_idx+1:row_idx+d, end-d+1:end) = eye(d);
            row_idx = row_idx + d;
        end
        col_idx = col_idx + k * d;
    end

    L = M'*M;
    L(end-d+1:end,1:end-d) = 0;
%     L = L\M';
    
    L = inv(L);
    LM = M';

    for i=1:numel(E)
        % Monomial basis of DOF for the i-th shape.
        x = x0(:,E{i});
        Qi = monomial_basis(x, x0_com, order);
        
        % Build 'B' matrices (Aqq in shape matching)
        [U, S, V] = svd(Qi*Qi');
        S = diag(S);
        S = 1./ max(S, 1e-4);
        Bi = Qi' * (V * diag(S) * U');
        B{i} = Bi;
    end
end

