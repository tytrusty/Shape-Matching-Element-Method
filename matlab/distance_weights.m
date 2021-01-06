function w = distance_weights(parts, X, alpha, beta)
    m=size(X,1);
    n=numel(parts);

    if nargin < 4
        beta = 1; % specifies a minimum cutoff distance 
    end
    
    % Weights
    w = zeros(m,n);

    % Compute per-shape distance weights
    for i = 1:n
        if isfield(parts{i}, 'hires_T')
            FV.faces = parts{i}.hires_T;
            FV.vertices = parts{i}.hires_x0';
        else
            FV.faces = parts{i}.T;
            FV.vertices = parts{i}.x0';
        end

        [dist, surf_X] = point2trimesh(FV, 'QueryPoints', X, 'UseSubSurface', false);
        
        rays = X - surf_X;
        rays = rays ./ vecnorm(rays,2,2);
        
        % Shoot rays from each surface point through the query points.
        % If the hit points are within the cutoff distance, use the
        % distance between the surface point (surf_X) and the hit point.
        min_ray_dist = inf(size(X,1),1);
        for j = [ [1:i-1] [i+1:n] ]
            if isfield(parts{j}, 'hires_T')
                faces = parts{j}.hires_T;
                verts = parts{j}.hires_x0';
            else
                faces = parts{j}.T;
                verts = parts{j}.x0';
            end
            P0=verts(faces(:,1),:);
            P1=verts(faces(:,2),:);
            P2=verts(faces(:,3),:);
            
            for k = 1:size(X,1)
                or = surf_X(k,:);
                D = rays(k,:);
                [dist_ray,~] = ray_tri(P0, P1, P2, or, D);
                if numel(dist_ray) > 0
                    dist_ray = min(dist_ray);
                    min_ray_dist(k) = min(dist_ray, min_ray_dist(k));
%                     if k == 55
%                    plot3(X(k,1),X(k,2),X(k,3),'.','Color','g','MarkerSize',20);
%                    plot3(surf_X(k,1),surf_X(k,2),surf_X(k,3),'.','Color','b','MarkerSize',30);
%                    hitpnt = or + dist_ray * D;
%                    plot3(hitpnt(1),hitpnt(2),hitpnt(3),'o','Color','r','MarkerSize',30);
%                    plot3(P0(fids(1),1),P0(fids(1),2),P0(fids(1),3),'o','Color','r','MarkerSize',30);                
%                     end
                end
            end
        end
    
        bound = max(beta, min(alpha, min_ray_dist));
        w(:,i) = max(1 -  (abs(dist) ./ bound), 0);
    end

end