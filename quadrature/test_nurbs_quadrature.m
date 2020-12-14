function test_nurbs_quadrature
fig=figure(1);
parts=nurbs_from_iges('trident.iges',14,0);
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
[bins, face_bins, OT] = octree_mesh(faces, verts, 20);

%%%%%%%%%%%%%%%%%%%%%%%
% Drawing each triangle surf
for i=1:numel(face_bins)
%     trisurf(face_bins{i}, bins{i}(:,1),bins{i}(:,2),bins{i}(:,3), ...
%             'FaceAlpha',0.4,'EdgeColor','none');
    trisurf(face_bins{i}, verts(:,1),verts(:,2),verts(:,3), ...
    'FaceAlpha',0.4,'EdgeColor','none');
    hold on;
end
xlabel('x'); ylabel('y'); zlabel('z');
axis equal

% Example intersection ray
y= 5; z = 14.5;
valid_bins = find(all(([y z] > OT.BinBoundaries(:,2:3) & ...
                       [y z] < OT.BinBoundaries(:,5:6)),2));
      
% Plotting intersection bins             
boxH = OT.plot('Visible',0);
for ii = 1:numel(valid_bins)
   i = valid_bins(ii);
   set(boxH(i),'Color','r','LineWidth', 1,'Visible',1)
   hold on;
end

xval=40;
plot3([-xval; xval], [y; y], [z; z],'Color','b', 'LineWidth',3);
hold on;
axis image, view(3)

face_intersect = face_bins(valid_bins);
face_intersect = face_intersect(~cellfun(@isempty, face_intersect));

vert_intersect = bins(valid_bins);
vert_intersect = vert_intersect(~cellfun(@isempty, vert_intersect));

% fids=[];
%intersect_faces 
for ii = 1:numel(vert_intersect)
    vert_ii = vert_intersect{ii};
    face_ii = face_intersect{ii};
%     P0=vert_ii(face_ii(:,1),:);
%     P1=vert_ii(face_ii(:,2),:);
%     P2=vert_ii(face_ii(:,3),:);
    P0=verts(face_ii(:,1),:);
    P1=verts(face_ii(:,2),:);
    P2=verts(face_ii(:,3),:);
    or=[-xval y z];
    D = [xval y z] - or;
    D = D ./ norm(D);
    [dist,fids] = ray_tri(P0, P1, P2, or, D);
    %trisurf(face_ii(fids,:), vert_ii(:,1),vert_ii(:,2),vert_ii(:,3),'FaceColor','b','FaceAlpha',1.0);
    trisurf(face_ii(fids,:), verts(:,1),verts(:,2),verts(:,3),'FaceColor','b','FaceAlpha',1.0);
%     trisurf(face_ii, vert_ii(:,1),vert_ii(:,2),vert_ii(:,3),'FaceColor','b','FaceAlpha',1.0);
    hold on;
end

end