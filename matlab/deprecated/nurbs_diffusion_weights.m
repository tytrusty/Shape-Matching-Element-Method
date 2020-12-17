function [a] = nurbs_diffusion_weights(parts, X)

    function xi = GetX(Ji,qi)
        Ji = reshape(Ji,[],3,size(Ji,4));
        Ji = permute(Ji,[3 2 1]);
        xi = squeeze(sum(Ji .* qi,1));
    end

    m=size(X,2);
    n=numel(parts);
    
    % Weights
    a = zeros(m,n);
    res=20;
    N=res*res;
    
    % Solving diffusion for each edge
    for i = 1:n
%         u = linspace(parts{i}.u_range(1),parts{i}.u_range(2),res)';
%         v = linspace(parts{i}.v_range(1),parts{i}.v_range(2),res)';
        % Gauss-Legendre shouldn't be used when the patch is period
        % (cylinder,sphere,etc.)
        [u,~]=lgwt(res, parts{i}.u_range(1)+1e-4, parts{i}.u_range(2)-1e-4);
        [v,~]=lgwt(res, parts{i}.v_range(1)+1e-4, parts{i}.v_range(2)-1e-4);
        [u_g,w_u]=lgwt(res, parts{i}.u_range(1), parts{i}.u_range(2));
        [v_g,w_v]=lgwt(res, parts{i}.v_range(1), parts{i}.v_range(2));
        [w_u,w_v] = meshgrid(w_u,w_v);
        w = w_u .* w_v;
        w = w(:);
        
        [q, obj_J] = nurbs_coords(parts{i}.nurbs, u, v);
        x = GetX(obj_J, q);
        
        [q, obj_J] = nurbs_coords(parts{i}.nurbs, u_g, v_g);
        x_g = GetX(obj_J, q);
        

        f = ones(N,1);
        L=zeros(N,N);
        inv_fourpi = 1 / (4 * pi);
        for j = 1:N
            diff = x(:,j) - x_g;

            % Single layer formulation
            L(j,:) = inv_fourpi * (1 ./ vecnorm(diff)) .* w';
        end
        mu = L \ f;
        mu = mu';
        
        u = zeros(m,1);
        for j = 1:m
            diff = X(:,j) - x_g;

            % Single layer
            green = inv_fourpi * (1 ./ vecnorm(diff)) .* mu .* w';
            u(j) = sum(green);
        end
        max(u)
        min(u)
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

