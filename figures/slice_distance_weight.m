clear;
iges_file = 'mug.iges';
parts = nurbs_from_iges(iges_file);
[V, F] = triangulate_iges(parts);

% define the boundary of the cross section to be a rectangle 
U_sqr = [min(V(:,1))-5 min(V(:,3))-5;...
         min(V(:,1))-5 max(V(:,3))+5;...
         max(V(:,1))+5 max(V(:,3))+5;...
         max(V(:,1))+5 min(V(:,3))-5];
E_sqr = [1 2; 2 3; 3 4; 4 1];
writePOLY_triangle('./test_sqr.poly', U_sqr, E_sqr, []);
[TV_sqr,TF_sqr,TN_sqr] = triangle('./test_sqr.poly');

alpha = 1.0; % cutoff distance

step = 50; % number of the slicing steps 
TF_list = cell(step,1);
X_list = cell(step,1);
w_list = cell(step,1);

ii = 1;
for d = linspace(min(V(:,2))+0.2,max(V(:,2)-0.2),step)
  plane = [0 1 0 -d]; % define the plane
  [U,E,J] = slice_triangles(V,F,plane); % find the cross section of the nurbs
  U_seg = [U(:,1) U(:,3)];
  
  % concat the rectange and the cross section contour
  U = [U_seg; U_sqr];
  E = [E; E_sqr + repmat(max(E,[],'all'),size(E_sqr,1),size(E_sqr,2))];
  writePOLY_triangle('./test.poly', U, E, []);
  
  % triangulate the cross section
  [TV,TF,TN] = triangle('./test.poly','Quality',true,'MaxArea',0.1);
    
  % compute the distance weight
  X = [TV(:,1) repmat(d,size(TV,1),1) TV(:,2)];
  w = distance_weights(parts, X, alpha, false);
  
  TF_list{ii} = TF;
  X_list{ii} = X;
  w_list{ii} = w;
  
  ii = ii + 1
end

% for visualization
clf;
figure(1);
t = tsurf(F,V,'FaceColor',[0.8 0.8 0.8],'FaceAlpha',0.4,'EdgeColor','none');
hold on;
s = tsurf([1 1 1],V,'EdgeColor','none');
colormap summer
colorbar
caxis([0 1])
axis equal
hold off;
for ii = 1:step
  set(s, 'Faces', TF_list{ii}, 'Vertices', X_list{ii}, 'CData', w_list{ii}(:,2), 'FaceAlpha', 0.8);
  drawnow;
end

% triangulate the nurbs parts
function [verts,faces] = triangulate_iges(parts)
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
  % check and remove degenerate faces
  dblA = doublearea(verts,faces);
  deg_fid = find(dblA <= 0);
  faces(deg_fid, :) = [];
end