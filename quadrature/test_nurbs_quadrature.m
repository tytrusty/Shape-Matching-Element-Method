function test_nurbs_quadrature
fig=figure(1);
clf;
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

%clf;
%trisurf(faces, verts(:,1),verts(:,2),verts(:,3));
clf;
OT = OcTree(verts);
% OT = OcTree(verts,'binCapacity',10);
%        boxH = OT.plot;
%        cols = lines(OT.BinCount);
%        doplot3 = @(p,varargin)plot3(p(:,1),p(:,2),p(:,3),varargin{:});
%        for i = 1:OT.BinCount
%            set(boxH(i),'Color',cols(i,:),'LineWidth', 1+OT.BinDepths(i))
%            doplot3(verts(OT.PointBins==i,:),'.','Color',cols(i,:))
%        end
%        axis image, view(3)

% Column vector of vertex indices 1:N
ids = (1:size(verts,1))';

% assign points to bins then collapse
bins = cell(OT.BinCount,1);
bin_ind = cell(OT.BinCount,1);
for i = 1:numel(bins)
    bins{i} = verts(OT.PointBins==i,:);
    bin_ind{i} = ids(OT.PointBins==i,:);
end
% Place faces into octree bins
face_bins = cell(size(bins,1),1);
newIDBins = cell(size(bins,1),1);
nfaces = size(faces,1);
all_faces = vertcat(faces(:,1), faces(:,2), faces(:,3));
for i=1:numel(bins)
    logi = ismember(all_faces, bin_ind{i});
    logicut = logi(1:nfaces) + logi(nfaces+1:nfaces*2) + logi((nfaces*2)+1:end);
    logicut(logicut > 0)=1;
    
    % Adds the vertex set of each face that intersects with this bin.
    face_bins{i} = faces(logical(logicut),:); 
    
    % For all faces that interesct with this bin, add all of their
    % associated vertices + vertex ids to the vertex bin and vertex id bin
    inpoints = unique(face_bins{i});
    inpoints(ismember(inpoints, bin_ind{i})) = [];     % must not be in set
    bins{i} = vertcat(bins{i}, verts(inpoints,:));     % add new vertices
    bin_ind{i} = vertcat(bin_ind{i}, ids(inpoints,:)); % add new vert ids
    
    newIDBins{i} = (1:size(bin_ind{i},1))';
    bin_size = size(face_bins{i},1);
    face_stack = vertcat(face_bins{i}(:,1), face_bins{i}(:,2), face_bins{i}(:,3));
    
    % Remap global vertex ids to local ids within the single bin.
    % These IDS in newIDBins serve as indices into bin_ind, which tell
    % us the global vertex ids.
    tempInds = arrayfun(@(x) find(bin_ind{i} == x,1,'first'), face_stack);
    face_stack = newIDBins{i}(tempInds);
    
    % face_bins now contains local indices.
    face_bins{i} = [face_stack(1:bin_size)  face_stack(bin_size+1:bin_size*2)  face_stack((bin_size*2)+1:end)];
end

%%%%%%%%%%%%%%%%%%%%%%%
% Drawing each triangle surf
for i=1:numel(face_bins)
    trisurf(face_bins{i}, bins{i}(:,1),bins{i}(:,2),bins{i}(:,3),'FaceAlpha',0.4,'EdgeColor','none');
    hold on;
end
xlabel('x'); ylabel('y'); zlabel('z');
axis equal
%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%
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
    P0=vert_ii(face_ii(:,1),:);
    P1=vert_ii(face_ii(:,2),:);
    P2=vert_ii(face_ii(:,3),:);
    or=[-xval y z];
    D = [xval y z] - or;
    D = D ./ norm(D);
    [dist,fids] = ray_tri(P0, P1, P2, or, D);
    trisurf(face_ii(fids,:), vert_ii(:,1),vert_ii(:,2),vert_ii(:,3),'FaceColor','b','FaceAlpha',1.0);
%     trisurf(face_ii, vert_ii(:,1),vert_ii(:,2),vert_ii(:,3),'FaceColor','b','FaceAlpha',1.0);
    hold on;
end

end