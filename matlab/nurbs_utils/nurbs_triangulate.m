function [T,UV] = nurbs_triangulate(nurbs, use_triangle)
    if nargin < 3
        use_triangle = 0;
    end

    if nurbs.is_trimmed
        if use_triangle
            [T,UV] = trimmed_trimesh_triangle(nurbs.line_1, ...
                nurbs.line_2, nurbs.line_N);
        else
            warning('Triangulating trimmed nurbs without triangle!');
            T = trimmed_trimesh(nurbs.line_1, nurbs.line_2, ...
                                nurbs.line_N, nurbs.UV);
            UV = nurbs.UV;
        end
        
    else
        % TODO -- make untrimmed_res optional so that you can optionally
        % add more untrimmed samples
        T = delaunay(nurbs.UV(1,:), nurbs.UV(2,:));
        UV = nurbs.UV;
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

    function [T,V] = trimmed_trimesh_triangle(p1, p2, N)
        offset = 0.2;	% offset along normal to generate hole points.
        midpoints = (p1+p2) ./ 2;               % line midpoints
        p_outside = midpoints + offset .* N;    % points outside surface

        % Define a boundary around the surface. This is just some arbitrary
        % number high enough so that the UV boundary will enclose all the
        % points. This UV boundary lets us define the points that are
        % outside the surface as holes for triangle to indicate where not
        % to triangulate.
        uv_bnd = 50;

        % Planar Straight Line Graph (PSLG) vertices for triangle.
        TUV = [p1 p2];
        nlines = size(p1,2);

        % Basic heuristic so that triangle density is reasonable.
        max_area_divisor = 100;
        max_area = min(range(TUV, 2)) / max_area_divisor;

        % Beginning index for UV boudnary points.
        bnd_idx = size(TUV,2) + 1;

        % Adding UV boundary vertices.
        TUV = [TUV [-uv_bnd -uv_bnd]' [uv_bnd -uv_bnd]' ...
                     [uv_bnd uv_bnd]' [-uv_bnd uv_bnd]'];

        % PSLG segments for triangle.
        E = [1:nlines; nlines+1:2*nlines]';

        % Add edges on the UV boundary.
        E = [E; [bnd_idx   bnd_idx+1; ...
                 bnd_idx+1 bnd_idx+2; ...
                 bnd_idx+2 bnd_idx+3; ...
                 bnd_idx+3 bnd_idx]];

        H = p_outside';

        % Invoke triangle. This assumes your path is set is set in
        % path_to_triangle file in GPToolbox in the file:
        % gptoolbox\wrappers\path_to_triangle.m
        %
        % -j Flags removes vertices not included in the triangulation.
        [V,T] = triangle(TUV', E, H, 'Quality','MaxArea', max_area, ...
            'Flags','-j');
        V = V';
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