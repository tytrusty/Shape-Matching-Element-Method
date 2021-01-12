function test_nurbs_quadrature
fig=figure(1);

parts=nurbs_from_iges('rounded_cube.iges',8,0);
parts=nurbs_plot(parts);

X = raycast_quadrature(parts, [8 8], 5);

%[flag, t, lambda] = ray_mesh_intersect(ray(:,1:3), ray(:,4:6), fv.vertices, fv.faces);

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

clf;

% Create spatial structure
[face_bins, OT] = octree_mesh(faces, verts);

%%%%%%%%%%%%%%%%%%%%%%%
% Drawing each triangle surf
for i=1:numel(face_bins)
    trisurf(face_bins{i}, verts(:,1),verts(:,2),verts(:,3), ...
    'FaceAlpha',0.6,'EdgeColor','none');
    hold on;
end
xlabel('x'); ylabel('y'); zlabel('z');
axis equal

xrange = [OT.BinBoundaries(1,1) OT.BinBoundaries(1,4)];
yrange = [OT.BinBoundaries(1,2) OT.BinBoundaries(1,5)];
zrange = [OT.BinBoundaries(1,3) OT.BinBoundaries(1,6)];

subd = [8 8];
samples = 5;
[U, V] = meshgrid(linspace(yrange(1), yrange(2), subd(1)), ...
                  linspace(zrange(1), zrange(2), subd(2)));
UV = [U(:) V(:)];

%     or = [xrange(2) coord];
%     
%     D = [xrange(1) coord] - or;
%     D = D ./ norm(D);
%     
%   
%[flag, t, lambda] = ray_mesh_intersect(ray(:,1:3), ray(:,4:6), fv.vertices, fv.faces);
ray_origins = [repelem(xrange(2), size(UV,1), 1) UV];
ray_dirs = repmat([-1 0 0], size(UV,1), 1);

% TODO this wont suffice since we need all intersections :/
% [flag, t, lambda] = ray_mesh_intersect(ray_origins, ray_dirs, verts, faces);

for i = 1:size(UV,1)
    % Intersect along x-axis
    [dist, face_intersect] = ray_intersect(xrange, UV(i,:), verts, face_bins, OT);
    trisurf(face_intersect, verts(:,1),verts(:,2),verts(:,3), ...
            'FaceColor','b', 'FaceAlpha', 1.0);
    hold on;
    
    if ~isempty(dist)
     	plot3(xrange, [UV(i,1); UV(i,1)], [UV(i,2); UV(i,2)],'Color','g', 'LineWidth',1);
        hold on;
    
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
                plot3(pnts(:,1),pnts(:,2),pnts(:,3),'.','Color','r', 'MarkerSize', 25)
                hold on;
            end
        else
            UV_rep = repelem(UV(i,:),numel(dist),1);
            pnts = xrange(2) - dist;
            pnts = [pnts UV_rep];
            plot3(pnts(:,1),pnts(:,2),pnts(:,3),'.','Color','r', 'MarkerSize', 15)
            hold on;
        end
        
    end
end
plot3(X(:,1),X(:,2),X(:,3),'o','Color','b', 'MarkerSize', 15)

end