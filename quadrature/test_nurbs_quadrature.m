function test_nurbs_quadrature
fig=figure(1);
parts=nurbs_from_iges('rocket_4.iges',14,0);
parts=plot_nurbs(parts);

faces=[];
verts=[];
for ii=1:numel(parts)
    fvc = surf2patch(parts{ii}.plt,'triangles');
    fvc.faces = fvc.faces + size(verts,1);
    faces=[faces; fvc.faces];
    verts=[verts; fvc.vertices];   
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

subd = [8 12];
samples = 5;
[U, V] = meshgrid(linspace(yrange(1), yrange(2), subd(1)), ...
                  linspace(zrange(1), zrange(2), subd(2)));
UV = [U(:) V(:)];
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
            plot3(pnts(:,1),pnts(:,2),pnts(:,3),'.','Color','r', 'MarkerSize', 25)
            hold on;
        end
        
    end
end

end