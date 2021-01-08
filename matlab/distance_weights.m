function w = distance_weights(parts, X, alpha, enable_secondary_rays)
    m=size(X,1);
    n=numel(parts);

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

        [dist, surf_X] = point2trimesh(FV, 'QueryPoints', X, ...
                                       'UseSubSurface', false);
        
        if enable_secondary_rays
            rays = X - surf_X;
            rays = rays ./ vecnorm(rays,2,2);

            % Shoot rays from each surface point through the query points.
            % If the hit points are within the cutoff distance, use the
            % distance between the surface point (surf_X) and the hit
            % point.
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
                    end
                end
            end
            bound = min(alpha, min_ray_dist);
        else
            bound = alpha;
        end
        w(:,i) = min(max(1 -  (abs(dist) ./ bound), 0), 1);
    end

end