function [W,W_I,W_S] = build_weight_matrix(w, d, k, varargin)
    p = inputParser;
    addParameter(p, 'Truncate', true);
    addParameter(p, 'Threshold', 1e-2);
    parse(p,varargin{:});
    
    truncate = p.Results.Truncate;
    threshold = p.Results.Threshold;
    
    m = size(w,1);  % number of points
    n = size(w,2);  % number of shapes
    
    W = cell(m,1);
    W_I = cell(m,1);
    W_S = cell(m,1);
    
    for i = 1:m
        
        if truncate
            ToKeep = w(i,:) > threshold;
        else
            ToKeep = ones(1,n);
        end
        
        % Shape indices for weights above the threshold.
        W_I{i} = find(ToKeep)';
        t = numel(W_I{i});
       
        % Forming selection and weight matrices.
        W{i} = zeros(d*k, d*k*t);
        W_S{i} = zeros(d*k*t, d*(k*n + 1));
        for j = 1:t
            I = W_I{i}(j);
            
            col_range = d*k*(j-1)+1:d*k*j; 
            W{i}(:, col_range) = eye(d*k)*w(i,I);
            
            
            row_range = d*k*(j-1)+1 : d*k*j;
            col_range = d*k*(I-1)+1:d*k*I; 
            W_S{i}(row_range, col_range) = eye(d*k);
        end     
    end
end