function [a] = compute_diffusion_weights(V, E, X)
    m=size(X,2);
    n=numel(E);
    
    % Weights
    a = zeros(m,n);
    
    % Solving diffusion for each edge
    for i = 1:n
        N = 8;
        % Generate points along line
        y = linspace(0,1, N);
        [y_g,w_y]=lgwt(N, 0,1);
        dir = V(:,E{i}(1)) - V(:,E{i}(2));
        y_g  = V(:,E{i}(2)) + dir * y_g';
        y  = V(:,E{i}(2)) + dir * y;
        
        f = ones(N,1);
        L=zeros(N,N);
        inv_twopi = 1 / (2 * pi);
        for j = 1:N
            diff = y(:,j) - y_g;

            % Single layer formulation
            L(j,:) = -inv_twopi * log(vecnorm(diff)) .* w_y';
        end
        mu = L \ f;
        mu = mu';
        
        u = zeros(m,1);
        for j = 1:m
            diff = X(:,j) - y_g;

            % Single layer
            green = -inv_twopi * log(vecnorm(diff)) .* mu .* w_y';
            u(j) = sum(green);
        end
        u = rescale(u);
%         u = max(u,0);
        % cmap=colormap(hot(256));
        % C = u ./ max(u);
        % u = max(u,0);
        %C = floor(rescale(u)*255) + 1;
        %C = cmap(C,:);
        %scatter(X(1,:),X(2,:),20,C,'filled');
        a(:,i) = u;
    end

    % Compute weighting coefficients
    options = optimoptions('quadprog','Display', 'off');

    for i = 1:m
        % The nearest points on each shape to our query point, d
        H = diag(1./a(i,:));
        H(isinf(H)) = 1e8;
        f = zeros(n,1);
        
        Aeq = ones(1,n);
        beq = 1;
        
        ub = ones(n,1);
        lb = zeros(n,1);
        a(i,:) = quadprog(H, f, [], [], Aeq, beq, lb, ub, [], options);
    end    
end

