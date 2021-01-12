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

% Create spatial structure
[face_bins, OT] = octree_mesh(faces, verts);
xrange = [OT.BinBoundaries(1,1) OT.BinBoundaries(1,4)];
yrange = [OT.BinBoundaries(1,2) OT.BinBoundaries(1,5)];
zrange = [OT.BinBoundaries(1,3) OT.BinBoundaries(1,6)];

% TODO get rid of OT stuff :)

half_inv_du = 1 / ray_samples(1) / 2;
half_inv_dv = 1 / ray_samples(2) / 2;
[U, V] = meshgrid(linspace(yrange(1), yrange(2), 1+ray_samples(1)), ...
                  linspace(zrange(1), zrange(2), 1+ray_samples(2)));
U = U(1:end-1, 1:end-1);
V = V(1:end-1, 1:end-1);
UV = [U(:)+half_inv_du V(:)+half_inv_dv];

dA = ((yrange(2)-yrange(1))/ray_samples(1)) * ...
     ((zrange(2)-zrange(1))/ray_samples(2)); 
 
X = [];
vol = [];

x_origin = xrange(2) + 1;
origins = [repelem(x_origin, size(UV,1), 1) UV];
ray_dirs = repmat([-1 0 0], size(UV,1), 1);

[t, nhits] = ray_mesh_intersections(origins, ray_dirs, verts, faces);
hit_ids = find(nhits > 0);

for i = 1:numel(hit_ids)
    idx = hit_ids(i);
    dist = t(idx,1:nhits(idx));
    
    % If distance list is even, then we split intersections into pairs
    % and sample between each intersection. Otherwise, we only create
    % a single point sample at each intersection.
    if rem(numel(dist),2) == 0
        dist = reshape(dist,2,[]);
        UV_rep = repelem(UV(idx,:),samples,1);
        for j = 1:size(dist,2)
            x_samples = linspace(dist(1,j),dist(2,j),samples+1);
            x_samples = x_samples(1:end-1);

            % Length of each segment sampled on the ray.
            dt = (dist(2,j)-dist(1,j)) / samples;

            pnts = x_origin - x_samples - dt/2;
            pnts = [pnts' UV_rep];
            vol_j = dt * dA;
            X = [X; pnts]; % inefficient, yea
            vol = [vol; repelem(vol_j, samples,1)];
        end
    else
        UV_rep = repelem(UV(idx,:),numel(dist),1);
        pnts = x_origin - dist;
        pnts = [pnts' UV_rep];
        X = [X; pnts];
        vol = [vol; repelem(dA, numel(dist),1)];
    end
end
X = X';

end