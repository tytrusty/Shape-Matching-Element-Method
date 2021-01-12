function  test_trimmed_nurbs_2

% Okay, V1 sucked! Simply discarding points is poor idea, so either need
% to reparameterize UV using something like least squares conformal map
% or do a better job of sampling on the trimmed surface. The latter is
% easier, so I'm gonna do that here with the ray casting integration
% strategy (same one I use for doin quadrature over volumes)
clf;
% data=iges2matlab('cylinder2.iges');
data=iges2matlab('wrench.iges');
% data=iges2matlab('wrench.igs');


% For hires mesh do the following:
%   - add p1 (plus last element of p2) to UV
%   - increase nlines in boundary_lines for line type
%   - in sample_ray make sure to exclude endpoints
%   - also kill degenerate triangles by checking area


function dist = point_line_min_distances(origins, p1, p2)
    % ref: http://paulbourke.net/geometry/pointlineplane/
    diff_X = origins(1,:)' - p1(1,:);
    diff_Y = origins(2,:)' - p1(2,:);
    dp = p2 - p1;
    norm_sqr = dot(dp,dp,1);

    u = (diff_X .* dp(1,:) + diff_Y .* dp(2,:)) ./ norm_sqr;
    u(:) = max(min(u(:), 1), 0);
    
    min_X = p1(1,:) + u .* dp(1,:);
    min_Y = p1(2,:) + u .* dp(2,:);
    
    diff_X = origins(1,:)' - min_X;
    diff_Y = origins(2,:)' - min_Y;
    dist = diff_X.^2 + diff_Y.^2;
end

function [p1,p2] = boundary_lines(boundary_srf)
    p1 = [];
    p2 = [];
    
    for i=1:boundary_srf.n
        boundary = data{boundary_srf.bdpt(i)};
        for j=1:numel(boundary.pscpt)
            curve = data{boundary.pscpt{j}};
            if curve.type == 110 % Line
                nlines =1;
            elseif curve.type == 126 % NURBS CURVE
                nlines = 50;
            end
            
            tt = linspace(curve.v(1), curve.v(2), nlines+1);
            p = nrbeval(curve.nurbs, tt); 
            p1 = [p1 p(1:2,1:end-1)];
            p2 = [p2 p(1:2,2:end)];         
        end
    end
    plot([p1(1,:); p2(1,:)],[p1(2,:); p2(2,:)], 'LineWidth',5, 'Color',[0.5 0.5 0.5]);
    hold on;
end

function data=init_untrimmed(data)
    for i=1:numel(data)  
       if data{i}.type == 128
           data{i}.is_trimmed = 0;
       end
    end
end

function UV = sample_ray(p1, p2, v_vals, u_range, nrays, samples_per_ray)
    UV = [];
    for jj = 1:nrays
        % Find intersection points.
        isect = (p1(2,:) <= v_vals(jj) & p2(2,:) >= v_vals(jj)) ...
              | (p1(2,:) >= v_vals(jj) & p2(2,:) <= v_vals(jj));

        if nnz(isect) > 0
            origin = [u_range(1) v_vals(jj)]';
            p1_isect = p1(:,isect);
            p2_isect = p2(:,isect);

            % reference for ray-line intersection:
            % https://rootllama.wordpress.com/2014/06/20/ray-line-segment-intersection-test-in-2d/
            v1 = origin - p1_isect;
            v2 = p2_isect - p1_isect;
            v3 = repmat([0 1]', 1, size(v1,2));

            t = (v2(1,:).*v1(2,:) - v2(2,:).*v1(1,:)) ./ dot(v2,v3,1);
            t = sort(t);
            t = uniquetol(t);
            assert(rem(numel(t),2) == 0, ...
                'Number of intersection is not even! Trimming failed.')

            t = reshape(t,2,[]);

            t_pnts = [];
            for kk = 1:size(t,2)
                t_shift = [t(1,kk)+1e-2 t(2,kk)-1e-2];
                if t_shift(1) < t_shift(2)
                    t_samples = linspace(t_shift(1),t_shift(2), samples_per_ray);
                    % Exclude endpoints
                    t_samples = t_samples(1,2:end-1);
                    new_t_pnts = origin + t_samples .* [1 0]';
                    t_pnts = [t_pnts new_t_pnts];
                end
            end

%             plot(p1_isect(1,:), p1_isect(2,:), '.', 'Color', 'm', 'MarkerSize', 10);
%             plot([u_range(1) u_range(2)], [v_vals(jj) v_vals(jj)],'-','Color','b','LineWidth',2);
%             plot([p1(1,isect); p2(1,isect)], [p1(2,isect); p2(2,isect)], '-', 'Color','r','LineWidth',2);
%             if numel(t_pnts) > 0
%                 plot(t_pnts(1,:), t_pnts(2,:), '.', 'Color', 'g', 'MarkerSize', 20);
%             end
            UV = [UV t_pnts];
        end
    end
end

function [UV,T] = trimmed_trimesh(p1, p2, v_vals, u_range, nrays, samples_per_ray)
    UV = sample_ray(p1, p2, v_vals, u_range, nrays, samples_per_ray);

    T = delaunay(UV(1,:), UV(2,:));
    
    z = zeros(1, size(T,1));
    Tu = UV(:,T(:,2)) - UV(:,T(:,1));
    Tv = UV(:,T(:,3)) - UV(:,T(:,1));
    T_double_area = cross([Tu; z], [Tv; z]);
    T_double_area = T_double_area(3,:);
    T(T_double_area < 1e-4, :) = [];
    

    % Purge triangles outside of the patch.
    centroids = (UV(:,T(:,1)) + UV(:,T(:,2)) + UV(:,T(:,3))) / 3;
    % plot(centroids(1,:), centroids(2,:), 'x', 'Color','r', 'MarkerSize',10);

    % Compute normals for each line segment sampled on boundary
    line = p2-p1;
    zero_z = repelem(0, size(p1,2));
    unit_z = repmat([0 0 1]',1, size(p1,2));
    N = cross([line; zero_z],unit_z);
    
    N = N ./ vecnorm(N);
    N=N(1:2,:);
        
%     Visualize normals
    pc = (p1+p2) ./ 2;
    p3 = pc + 0.2 .* N(1:2,:);
%     plot([pc(1,:); p3(1,:)], [pc(2,:); p3(2,:)],'b');

    % Compute minimum distance point from centroid to boundary lines 
    dist = point_line_min_distances(centroids, p1, p2);
    [~,I] = min(dist, [], 2);
    
    %plot([centroids(1,:); pc(1, I)], [centroids(2,:); pc(2, I)],'r');

    % Check angles between nearest point ray and normal.
    diff_p = centroids - p1(:, I);
    angles = dot(diff_p, N(:,I));
    TtoRemove = angles > 0;
    T(TtoRemove,:) = [];
    
    TUV = [p1 p2];
    nlines = size(p1,2);
    max_area = min(range(TUV, 2)) / 100
    
    bnd_idx = size(TUV,2) + 1;
    bnd = 50;
    TUV = [TUV [-bnd -bnd]' [bnd -bnd]' [bnd bnd]' [-bnd bnd]']; 
    E = [1:nlines; nlines+1:2*nlines]';
    
    
    E = [E; [bnd_idx bnd_idx+1; bnd_idx+1 bnd_idx+2; bnd_idx+2 bnd_idx+3; bnd_idx+3 bnd_idx]];
    H = p3'; %centroids(:,TtoRemove)';
%     [TV1,TF2,~] = triangle(UV', 'Quiet',false);

    plot([p1(1,:); p2(1,:)],[p1(2,:); p2(2,:)], 'LineWidth',2, 'Color',[0.5 0.5 0.5]);
    hold on;
%     plot([ TUV(1,E(:,1)); TUV(1,E(:,2))], [TUV(2,E(:,1)); TUV(2,E(:,2))], 'LineWidth',2, 'Color', 'r');
    plot(H(:,1)', H(:,2)', '.', 'Color','r', 'MarkerSize',30);
    hold on;
    
    
    [TV,TF] = triangle(TUV',E, H , 'Quality','MaxArea', max_area, 'Quiet',true);
    patch('Faces',TF,'Vertices',TV,'FaceColor','g', 'FaceAlpha', 0.2);
    axis equal;
end

% function 
data=init_untrimmed(data);

for ii=1:numel(data)   
    if data{ii}.type == 143
        clf;
    	srf_boundary = data{ii};
        [p1,p2] = boundary_lines(srf_boundary);
        
        % Sample a bunch of points on the surface
        srf_ii = srf_boundary.sptr;
                    
        % Mark surface as trimmed
        data{srf_ii}.is_trimmed = 1;
        
        nrays = 10;
        samples_per_ray = 10;
                
        v_range = [min([p1(2,:) p2(2,:)]) max([p1(2,:) p2(2,:)])];
        v_vals = linspace(v_range(1)+1e-4, v_range(2)-1e-4, nrays);
        u_range = [data{srf_ii}.u(1) data{srf_ii}.u(2)];
        [UV, T] = trimmed_trimesh(p1, p2, v_vals, u_range, nrays, samples_per_ray);
        data{srf_ii}.T = T;
        data{srf_ii}.UV = UV;
        data{srf_ii}.trim_lines = {p1, p2};
        patch('Faces',T,'Vertices',UV','FaceColor','g', 'FaceAlpha', 0.2);
        disp('Press a button to look at the next thing');
        waitforbuttonpress;
    end
end

clf;
for ii=1:numel(data)   
    if data{ii}.type == 128 && ~data{ii}.is_trimmed
        data{ii}.is_trimmed
        srf = data{ii};
        srf_res = 10;
        u = linspace(srf.u(1), srf.u(2), srf_res);
        v = linspace(srf.v(1), srf.v(2), srf_res);
        [U,V] = meshgrid(u,v);
        UV = [U(:) V(:)]';
        T = delaunay(UV(1,:), UV(2,:));
        
        patch('Faces',T,'Vertices',UV','FaceColor','g', 'FaceAlpha', 0.2);

        disp('test');
    end
end

end

