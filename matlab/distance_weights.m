function w = distance_weights(parts, X, alpha, enable_secondary_rays)
    m=size(X,1);
    n=numel(parts);

    % Weights
    w = zeros(m,n);

    % Compute per-shape distance weights
    for i = 1:n
        if isfield(parts{i}, 'hires_T')
            F = parts{i}.hires_T;
            V = parts{i}.hires_x0';
        else
            F = parts{i}.T;
            V = parts{i}.x0';
        end

        % Find closest points on the surface for each point.
        [sqrD,~,surf_X] = point_mesh_squared_distance(X,V,F);
        dist = sqrt(sqrD);

        if enable_secondary_rays
            rays = X - surf_X;
            rays = rays ./ vecnorm(rays,2,2);

            faces=[];
            verts=[];
            for j = [ [1:i-1] [i+1:n] ]
                if isfield(parts{j}, 'hires_T')
                    F = parts{j}.hires_T;
                    V = parts{j}.hires_x0';
                else
                    F = parts{j}.T;
                    V = parts{j}.x0';
                end
                F = F + size(verts,1);
                faces=[faces; F];
                verts=[verts; V];
            end

            % Shoot rays from each surface point through the query points.
            % If the hit points are within the cutoff distance, use the
            % distance between the surface point (surf_X) and the hit
            % point.
            [~, t] = ray_mesh_intersect(surf_X, rays, verts, faces);
            min_ray_dist = t;

            bound = min(alpha, min_ray_dist);
        else
            bound = alpha;
        end
        w(:,i) = min(max(1 -  (abs(dist) ./ bound), 0), 1);
    end

end