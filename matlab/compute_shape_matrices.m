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

    if strcmp(mode, 'global')
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
    elseif startsWith(mode,'local')
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

            if strcmp(mode, 'local')
                ATA_inv = inv(Ai'*Ai);
            elseif strcmp(mode, 'local_pinv')
                if r < size(Ai,2)
                    fprintf('deficient for shape %d\n', i);
                    ATA_inv = pinv(Ai'*Ai);
                else
                    ATA_inv = inv(Ai'*Ai);
                end
            elseif strcmp(mode, 'local_svd_truncated')
                epsilon = 1e-12;
                [U, S, V] = svd(Ai'*Ai);
                S = diag(S);
                t = find(abs(S) < epsilon, 1, 'first') - 1; %truncation point

                if isempty(t) % full rank
                    t = numel(S);
                end

                S = 1./ S(1:t);
                ATA_inv = V(:,1:t) * diag(S) * U(:,1:t)';
            end
            Li = ATA_inv * Ai' * (Sbi - Si*Ti);
            L(col_range,:) = Li;
        end
    end
end

