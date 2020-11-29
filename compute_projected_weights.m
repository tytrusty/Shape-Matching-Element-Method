function [a] = compute_projected_weights(V, E, X)

    m=size(X,2);
    n=numel(E);

    % Compute weighting coefficients
    a = zeros(m, n);
    options = optimoptions('lsqlin','Display', 'off');

    for i = 1:m
        % The nearest points on each shape to our query point, d
        C = zeros(size(X,1),numel(E));
        for j = 1:numel(E)
           Xj = V(:,E{j});
           dist = vecnorm((Xj-X(:,i)).^2);
           [~,I] = min(dist);
           C(:,j) = Xj(:,I);
           %C(:,j) = mean(V(:,E{j}),2);
        end
        
        % The query point
        d = X(:,i);

        % All weights non-negative
        A = -eye(n);
        b = zeros(n,1);

        % Weights sum to one.
        Aeq = ones(1,n);
        beq = 1;

        a(i,:) = lsqlin(C,d,A,b,Aeq,beq, [],[],[], options);
    end
    
end

