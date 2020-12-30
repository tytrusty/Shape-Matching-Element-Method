function [T,UV] = nurbs_triangulate(nurbs, untrimmed_res)

    use_triangle = 0;
    if nurbs.is_trimmed
        if use_triangle
            T=[];
        else
            T = trimmed_trimesh(nurbs.line_1, nurbs.line_2, ...
                                nurbs.line_N, nurbs.UV);
        end
        
    else
        u = linspace(nurbs.srf.u(1), nurbs.srf.u(2), untrimmed_res);
        v = linspace(nurbs.srf.v(1), nurbs.srf.v(2), untrimmed_res);
        [U,V] = meshgrid(u,v);
        UV = [U(:) V(:)]';
        T = delaunay(UV(1,:), UV(2,:));
    end

    function T = trimmed_trimesh(p1, p2, N, UV)
        T = delaunay(UV(1,:), UV(2,:));

        % Purge triangles outside of the patch.
        centroids = (UV(:,T(:,1)) + UV(:,T(:,2)) + UV(:,T(:,3))) / 3;

        % Compute minimum distance point from centroid to boundary lines 
        dist = point_line_min_distances(centroids, p1, p2);
        [~,I] = min(dist, [], 2);

        % Check angles between nearest point ray and normal.
        diff_p = centroids - p1(:, I);
        angles = dot(diff_p, N(:,I));
        TtoRemove = angles > 0;
        T(TtoRemove,:) = [];
    end

    function dist = point_line_min_distances(origins, p1, p2)
        % ref: http://paulbourke.net/geometry/pointlineplane/
        diff_X = origins(1,:)' - p1(1,:);
        diff_Y = origins(2,:)' - p1(2,:);
        dp = p2 - p1;
        norm_sqr = dot(dp,dp,1);

        t = (diff_X .* dp(1,:) + diff_Y .* dp(2,:)) ./ norm_sqr;
        t(:) = max(min(t(:), 1), 0);

        min_X = p1(1,:) + t .* dp(1,:);
        min_Y = p1(2,:) + t .* dp(2,:);

        diff_X = origins(1,:)' - min_X;
        diff_Y = origins(2,:)' - min_Y;
        dist = diff_X.^2 + diff_Y.^2;
    end
    
end