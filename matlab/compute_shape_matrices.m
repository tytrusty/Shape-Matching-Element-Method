function L = compute_shape_matrices(x0, x0_com, E, order, mode)
    if nargin < 5
        mode = 'default';
    end
    fprintf('Computing shape matrices in mode: %s', mode);
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
        row_ranges{i} = A_rows+1:A_rows + d*numel(E{i});
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

    if strcmp(mode, 'default')
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
    end
    %TODO still trying out stuff on local solves
%     elseif startsWith(str,'local')
%     	for i=1:numel(E)            
%             % Forming local system matrix
%             Ai = zeros(d*numel(E{i}), d*(k+1));
%             row_idx = 0;
%             for j=1:numel(E{i})
%                 for n=1:d
%                     col_start = (n-1)*k + 1;
%                     col_end   = n*k;
%                     Ai(row_idx+n, col_start:col_end) = M{i}(:,j);
%                 end
%                 Ai(row_idx+1:row_idx+d, end-d+1:end) = eye(d); % Identity
%                 row_idx = row_idx + d;
%             end
%             
%             Ai = A(row_ranges{i},:);
%             
%             r = rank(Ai);
%             fprintf('Local A for shape %d: nrows=%d ncols=%d rank=%d\n', ...
%                 i, size(Ai,1), size(Ai,2), r);
%             
%             ATAi_inv = inv(Ai'*Ai);
%             Li = (Ai'*Ai) \ Ai';
%             inv1 = inv(Ai'*Ai);
%             inv2 = pinv(Ai'*Ai);
% end
end

