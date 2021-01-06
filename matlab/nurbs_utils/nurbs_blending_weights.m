function w = nurbs_blending_weights(parts, X, alpha)

m=size(X,1);
n=numel(parts);

% Compute per-shape distance weights
w = distance_weights(parts, X, alpha);

% Compute weighting coefficients with constraints
options = optimoptions('quadprog','Display', 'off');

for i = 1:m
    % The nearest points on each shape to our query point, d
    H = diag(1./w(i,:));
    H(isinf(H)) = 1e8;
    f = zeros(n,1);

    Aeq = ones(1,n);
    beq = 1;

    ub = ones(n,1);
    lb = zeros(n,1);
    w(i,:) = quadprog(H, f, [], [], Aeq, beq, lb, ub, [], options);
end    

end