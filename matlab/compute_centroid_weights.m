function [a] = compute_weights(V,E,X)

    m=size(X,2);
    n=size(E,1);

    % Compute weighting coefficients
    a = zeros(m, n);
    options = optimoptions('lsqlin','Display', 'off');

    for i = 1:m
        C = V;
        d = X(:,i);

        % All weights non-negative
        A = -eye(size(V,2));
        b = zeros(size(V,2),1);

        % Weights sum to one.
        Aeq = ones(1,n);
        beq = 1;

        a(i,:) = lsqlin(C,d,A,b,Aeq,beq, [],[],[], options);
    end
    
end

