function L = compute_shape_matrices(x0, x0_com, E, order, mode)
    if nargin < 5
        mode = 'global';
    end
    fprintf('Computing shape matrices in mode: %s \n', mode);
    d = size(x0,1); % dimension

    % Compute number of monomials
    k = basis_size(d, order);
    
    % Monomial basis vectors
    M = cell(numel(E),1);
    
    % Computing size of A matrix
    A_rows = 0;
    row_ranges = cell(numel(E),1);
    for i=1:numel(E)
        x = x0(:,E{i});
        M{i} = monomial_basis(x, x0_com, order);
        row_ranges{i} = d*A_rows+1:d*A_rows + d*numel(E{i});
        A_rows = A_rows + numel(E{i});
    end
    A_cols = d*k*numel(E) + d;
    A = zeros(d * A_rows, A_cols);

    col_idx = 0;
    row_idx = 0;
    % Forming global system matrix ([A S] matrix) 
    for i=1:numel(E)
        for j=1:numel(E{i})
            for n=1:d
                col_start = col_idx + (n-1)*k + 1;
                col_end   = col_idx + n*k;
                A(row_idx+n, col_start:col_end) = M{i}(:,j);
            end
            A(row_idx+1:row_idx+d, end-d+1:end) = eye(d); % Identity term
            row_idx = row_idx + d;
        end
        col_idx = col_idx + k * d;
    end
    
    
    fprintf('A matrix: nrows=%d ncols=%d rank=%d\n', ...
        size(A,1), size(A,2), rank(A));

    ATA = A'*A;
    ATA(end-d+1:end,1:end-d) = 0; % applying center of mass constraint

    if strcmp(mode, 'hierarchical')
        n = size(x0,2);

        L = zeros(d*(k*numel(E)+1), numel(x0));

        % Setting constant term rows
        for i=1:d
            L(end-d+i,i:d:end) = 1/n;
        end

        I = eye(numel(x0));

        % 'S' matrix (stacked identities)
        S = repmat(eye(d),n,1);

        % Matrix that computes the mean of the nodal values.
        T = repmat(eye(d)*(1/n),1,n);

        k_linear = nchoosek(d+1-1,1);
        linear_cols = zeros(k_linear*d*numel(E),1);
        % Move this to function so it supports arbitrary orders
        for i = 1:numel(E)
            for j = 1:d
                range_b = k_linear*d*(i-1) + k_linear*(j-1) + 1;
                range_e = k_linear*d*(i-1) + k_linear*j;
                i_range = range_b:range_e;

                % (+1 should be + k_linear for order 2 or the cumulative
                % sum of the previous k values for higher orders
                col_start = k*d*(i-1)+ k*(j-1) + 1;
                col_end = k*d*(i-1)+ k*(j-1) + k_linear;
                col_range = col_start : col_end;
                linear_cols(i_range) = col_range;
            end
        end

        A_lin = A(:,linear_cols);
        L_lin = (A_lin'*A_lin) \ (A_lin' * (I - S*T));
        L(linear_cols,:) = L_lin;

        if order == 2
            k_quad = nchoosek(d+2-1,2);
            quad_cols = zeros(k_quad*numel(E),1);

            for i = 1:numel(E)
                for j = 1:d
                    range_b = k_quad*d*(i-1) + k_quad*(j-1) + 1;
                    range_e = k_quad*d*(i-1) + k_quad*j;
                    i_range = range_b:range_e;

                    col_start = k*d*(i-1)+ k*(j-1) + k_linear + 1;
                    col_end = k*d*(i-1)+ k*(j-1) + k_linear + k_quad;
                    col_range = col_start : col_end;
                    quad_cols(i_range) = col_range;
                end
            end

            A_quad = A(:,quad_cols);
            rhs = A_quad' * (I - S * T -  A_lin*L_lin);
            L_quad = (A_quad'*A_quad) \ rhs;
            L(quad_cols,:) = L_quad;
        end
    elseif strcmp(mode, 'global')
        L = ATA \ A';
        
    elseif strcmp(mode, 'global_pinv')
        L = pinv(ATA) * A';
        
    elseif strcmp(mode, 'global_svd_truncated')
        epsilon = 1e-12;
        [U, S, V] = svd(ATA);
        S = diag(S);
        t = find(abs(S) < epsilon, 1, 'first') - 1; %truncation point
        
        if isempty(t) % full rank
            t = numel(S);
        end
        
        S = 1./ S(1:t);
        truncated_inv = V(:,1:t) * diag(S) * U(:,1:t)';
        L = truncated_inv * A';

    elseif strcmp(mode,'local')
        n = size(x0,2);
        L = zeros(d*(k*numel(E)+1), numel(x0));

        % Setting constant term rows
        for i=1:d
            L(end-d+i,i:d:end) = 1/n;
        end

    	for i=1:numel(E)
            m = numel(E{i});

            % Selecting subset of b values for this shape.
            Sbi = zeros(d*m, numel(x0));
            Sbi(:, row_ranges{i}) = eye(d*m);

            % 'S' matrix (stacked identities)
            Si = repmat(eye(d),m,1);

            % Matrix that computes the mean of the nodal values.
            Ti = repmat(eye(d)*(1/n),1,size(x0,2));

            col_range = d*k*(i-1)+1:d*k*i;

            % Getting the local system for this shape.
            Ai = A(row_ranges{i},col_range);

            r = rank(Ai);
            fprintf('Local A for shape %d: nrows=%d ncols=%d rank=%d\n', ...
                i, size(Ai,1), size(Ai,2), r);

            % Inverting local block. If low rank use pseudoinverse.
            if r < size(Ai,2)
                fprintf('deficient for shape %d\n', i);
                ATA_inv = pinv(Ai'*Ai);
            else
                ATA_inv = inv(Ai'*Ai);
            end

            Li = ATA_inv * Ai' * (Sbi - Si*Ti);
            L(col_range,:) = Li;
        end
    end
end

