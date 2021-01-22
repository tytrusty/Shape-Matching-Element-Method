function [w,w_I] = nurbs_blending_weights(parts, X, alpha, varargin)
    p = inputParser;
    addParameter(p, 'Truncate', true);
    addParameter(p, 'Threshold', 1e-4);
    addParameter(p, 'Enable_Secondary_Rays', true);
    parse(p,varargin{:});

    truncate = p.Results.Truncate;
    threshold = p.Results.Threshold;
    enable_secondary_rays = p.Results.Enable_Secondary_Rays;

    m=size(X,1);
    n=numel(parts);
    w_I = cell(m,1);

    % Compute per-shape distance weights
    w = distance_weights(parts, X, alpha, enable_secondary_rays);

    % Compute weighting coefficients with constraints
    options = optimoptions('quadprog','Display', 'off');

    % Quadratic program for each point that solves for the final blending
    % weights.
    for i = 1:m
        H = diag(1./w(i,:));
        H(isinf(H)) = 1e8;
        f = zeros(n,1);

        Aeq = ones(1,n);
        beq = 1;

        ub = ones(n,1);
        lb = zeros(n,1);
        w(i,:) = quadprog(H, f, [], [], Aeq, beq, lb, ub, [], options);

        % If truncate is enabled, make all weights above the threshold.
        % This will simplify other terms, making things much faster.
        if truncate
            ToKeep = w(i,:) > threshold;
        else
            ToKeep = ones(1,n);
        end

        % Shape indices for weights above the threshold.
        w_I{i} = find(ToKeep)';
    end
end