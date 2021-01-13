function [X, vol] = raycast_quadrature(parts, ray_samples, samples)
% Note: assuming nurbs_plot has been called so that we can use the
% plots data in part.plt

% Triangulate nurbs patch
faces=[];
verts=[];
for ii=1:numel(parts)
    if isfield(parts{ii}, 'hires_T')
        F = parts{ii}.hires_T;
        V = parts{ii}.hires_x0';
    else
        F = parts{ii}.T;
        V = parts{ii}.x0';
    end
    F = F + size(verts,1);
    faces=[faces; F];
    verts=[verts; V];
end

% Finding the bounding box and storing the range in each direction.
min_bnd = min(verts,[],1);
max_bnd = max(verts,[],1);
xrange = [min_bnd(1) max_bnd(1)];
yrange = [min_bnd(2) max_bnd(2)];
zrange = [min_bnd(3) max_bnd(3)];

% Sampling ray origins along the YZ plane.
half_inv_du = (yrange(2)-yrange(1)) / ray_samples(1) / 2;
half_inv_dv = (zrange(2)-zrange(1)) / ray_samples(2) / 2;
[U, V] = meshgrid(linspace(yrange(1), yrange(2), 1+ray_samples(1)), ...
                  linspace(zrange(1), zrange(2), 1+ray_samples(2)));
U = U(1:end-1, 1:end-1);
V = V(1:end-1, 1:end-1);
UV = [U(:)+half_inv_du V(:)+half_inv_dv];

% Integration area around the ray.
dA = ((yrange(2)-yrange(1))/ray_samples(1)) * ...
     ((zrange(2)-zrange(1))/ray_samples(2)); 
 
% Generating and intersecting each ray against the model.
x_origin = xrange(2) + 1;
origins = [repelem(x_origin, size(UV,1), 1) UV];
ray_dirs = repmat([-1 0 0], size(UV,1), 1);

% Finding all hits along the rays.
[t, nhits, fids] = ray_mesh_intersections(origins, ray_dirs, verts, faces);
hit_ids = find(nhits > 0);

% Get X component of face normals. This is all we need to check in/out
% since our rays are X-axis aligned.
cross2 = @(a,b,c) a(:,2).*b(:,3)-a(:,3).*b(:,2);
NormalX = @(p1,p2,p3) cross2(p2-p1,p3-p1);

X = [];
vol = [];

% For each ray with 1 or more hits, sample along the intersection intervals
% inside the model. Using a CSG raytracing-like approach where we compute
% the union of intersection intervals (or at least try to :)).
for i = 1:numel(hit_ids)
    idx = hit_ids(i);
    dist = t(idx,1:nhits(idx));
    fids_i = fids(idx,1:nhits(idx));

    % Getting X component of the normal on each hit face.
    Nx = NormalX(verts(faces(fids_i,1),:), verts(faces(fids_i,2),:), ...
                 verts(faces(fids_i,3),:));
    in_rays = (Nx >= 0) + -1*(Nx < 0);
    intervals = {};
    
    % Construct intervals of intersection by unioning the objects. This
    % assumes things are closed, but is a fine solution for now.
    % If this works okay, I can quickly move it to c++.
    j=1;
    while j < numel(dist)+1
        % If first hit is an exit point, just take a point sample.
        if (in_rays(j) < 0)
            intervals{end+1} = dist(j);
        else
            inside_cnt = 1;
            t_begin = dist(j);
            while (inside_cnt ~= 0 && j < numel(dist))
                j = j+1;
                inside_cnt = inside_cnt + in_rays(j);
            end
            % Currently assuming each in ray is accompanied by an out ray.
            % There are plenty of cases where this will not happen, but
            % for now this will work.
            intervals{end+1} = [t_begin dist(j)];
        end
        j=j+1;
    end

    % Sample along each interval. If interval consists of only one point.
    % just create a single point sample (the volume will just be an area)
    for j=1:numel(intervals)
        if numel(intervals{j}) > 1
            UV_rep = repelem(UV(idx,:),samples,1);
            t_range = intervals{j};
            x_samples = linspace(t_range(1),t_range(2),samples+1);
            x_samples = x_samples(1:end-1);

            % Length of each segment sampled on the ray.
            dt = (t_range(2)-t_range(1)) / samples;

            pnts = x_origin - x_samples - dt/2;
            pnts = [pnts' UV_rep];
            vol_j = dt * dA;
            X = [X; pnts]; % inefficient, yea
            vol = [vol; repelem(vol_j, samples,1)];
        else
            pnts = x_origin - intervals{j};
            pnts = [pnts' UV(idx,:)];
            X = [X; pnts];
            vol = [vol; repelem(dA, numel(dist),1)];
        end
    end
end
X = X';

end