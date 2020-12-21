function [X, vol] = raycast_quadrature(parts, ray_samples, samples)
% Note: assuming nurbs_plot has been called so that we can use the
% plots data in part.plt

% Triangulate nurbs patch
faces=[];
verts=[];
for ii=1:numel(parts)
    fvc = surf2patch(parts{ii}.plt,'triangles');
    fvc.faces = fvc.faces + size(verts,1);
    faces=[faces; fvc.faces];
    verts=[verts; fvc.vertices];   
end

% Create spatial structure
[face_bins, OT] = octree_mesh(faces, verts);
xrange = [OT.BinBoundaries(1,1) OT.BinBoundaries(1,4)];
yrange = [OT.BinBoundaries(1,2) OT.BinBoundaries(1,5)];
zrange = [OT.BinBoundaries(1,3) OT.BinBoundaries(1,6)];


[U, V] = meshgrid(linspace(yrange(1), yrange(2), ray_samples(1)), ...
                  linspace(zrange(1), zrange(2), ray_samples(2)));
UV = [U(:) V(:)];      

dA = ((yrange(2)-yrange(1))/ray_samples(1)) * ...
     ((zrange(2)-zrange(1))/ray_samples(2));

X = [];
vol = [];
for i = 1:size(UV,1)
    % Intersect along x-axis
    [dist, ~] = ray_intersect(xrange, UV(i,:), verts, face_bins, OT);

    if ~isempty(dist)
        dist=sort(dist);
        
        % If distance list is even, then we split intersections into pairs
        % and sample between each intersection. Otherwise, we only create
        % a single point sample at each intersection.
        if rem(numel(dist),2) == 0
            dist = reshape(dist,2,[]);
            UV_rep = repelem(UV(i,:),samples,1);
            for j = 1:size(dist,2)
                x_samples = linspace(dist(1,j),dist(2,j),samples);
                pnts = xrange(2) - x_samples;
                pnts = [pnts' UV_rep];
                vol_j = (dist(2,j)-dist(1,j)) * dA;
                X = [X; pnts]; % inefficient, yea
                vol = [vol; repelem(vol_j, samples,1)];
            end
        else
            UV_rep = repelem(UV(i,:),numel(dist),1);
            pnts = xrange(2) - dist;
            pnts = [pnts UV_rep];
            X = [X; pnts];
            vol = [vol; dA];
        end
        
    end    
end
X = X';

end