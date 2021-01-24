function cutoff_distance=cutoff_heuristic(parts, T, n, X)
% Choose cutoff distance based on desired level of sparsity of blending
% weights.
%
% Given a target sparsity, T, heuristic searches for the distance
% cutoff that produces the desired level of sparsity, meaning that the
% ratio of nonzeros to total elements in the set of blending equals T.
%
% Arguments:
%   parts: cell array containing each part of model
%   T: the target sparsity ratio in (0, 1]
%   n (optional): the number of cutoff distance candidates to generate
%   X (optional): the set of points on which distance weights are computed
%
    if nargin < 3
        n = 20;
    end
    
    if nargin < 4
        YZ_samples = [9 9];
        X_samples = 5;
        [X, ~] = raycast_quadrature(parts, YZ_samples, X_samples);
        X = X';
    end

    % Find max possible cutoff
    x = [];
    for i=1:numel(parts)
       x = [x parts{i}.x0];
    end
    min_bnd = min(x, [], 2);
    max_bnd = max(x, [], 2);
    max_d = max(max_bnd-min_bnd);

    % candidate cutoff values
    cutoff_vals = linspace(0,max_d, n+1);
    cutoff_vals = cutoff_vals(2:end);
    
    % Binary search for clostest distance cutoff that satifies the
    % sparsity objective.
    L = 1;
    R = numel(cutoff_vals);
    eps = 1e-5;
    sparsity = zeros(numel(cutoff_vals),1);
    
    T = 1 - T;
    while L <= R
        index = floor((L + R) / 2);

        c = cutoff_vals(index);

        w = distance_weights(parts, X, c, true);
        sparsity(index) = sum(w(:) > eps) / numel(w);

        if sparsity(index) < T
            L = index + 1;
        elseif sparsity(index) > T
            R = index - 1;
        else
            index = m;
            break;
        end
    end
    cutoff_distance = cutoff_vals(index);
end