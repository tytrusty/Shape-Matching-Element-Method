clear;
iges_file = 'puft_no_collar.iges';
parts = nurbs_from_iges(iges_file);
[V, F, ~] = triangulate_iges(parts);

iges_file_wn = 'puft_no_collar.iges'; % for winding number
parts_wn = nurbs_from_iges(iges_file_wn);
[V_wn, F_wn] = triangulate_iges(parts_wn);

% define the boundary of the cross section plane to be a rectangle 
U_sqr = [min(V(:,1))-1 min(V(:,3))-1;...
         min(V(:,1))-1 max(V(:,3))+1;...
         max(V(:,1))+1 max(V(:,3))+1;...
         max(V(:,1))+1 min(V(:,3))-1];
E_sqr = [1 2; 2 3; 3 4; 4 1];
V_mean = mean(V);
hole_pos = [V_mean(:,1) V_mean(:,3)];
writePOLY_triangle('./test_sqr.poly', U_sqr, E_sqr, []);
[TV_sqr,TF_sqr,TN_sqr] = triangle('./test_sqr.poly');

tsurf(TF_sqr, TV_sqr);

alpha = 5.0; % cutoff distance

step = 3; % number of the slicing steps 
TF_list = cell(step,1);
X_list = cell(step,1);
w_list = cell(step,1);

ii = 1;
for d = linspace(min(V(:,2))+1.0,max(V(:,2)-1.0),step)
  plane = [0 1 0 -d]; % define the plane
  try
    [U,E,J] = slice_triangles(V,F,plane); % find the cross section of the nurbs
    U = [U(:,1) U(:,3)];
    U = [U; U_sqr];
    E = [E; E_sqr + repmat(max(E,[],'all'),size(E_sqr,1),size(E_sqr,2))];
  catch
    U = U_sqr;
    E = E_sqr;
  end
 
  % concat the rectange and the cross section contour
  writePOLY_triangle('./test.poly', U, E, []);
  
  tsurf(E, U);
  drawnow;
  
  % triangulate the cross section
  [TV,TF,TN] = triangle('./test.poly', 'Quality', true, 'MaxArea', 0.0001);
    
  % compute the distance weight
  X = [TV(:,1) repmat(d,size(TV,1),1) TV(:,2)];
  % remove the vertices outside the mesh
  W = winding_number(V_wn,F_wn,X);
  % visualize
  scatter3(X(:,1),X(:,2),X(:,3),10,(max(W,1e-15))./max(W));
  colorbar
  drawnow;
  
  w = distance_weights(parts, X, alpha, true);
  W = (max(W,1e-15))./max(W);
  w(W < 1e-1, :) = NaN;
  
  TF_list{ii} = TF;
  X_list{ii} = X;
  w_list{ii} = w;
  
  ii = ii + 1
end

%% for visualization
save_path = ['../output/slice/'];
mkdir(save_path);  
clf;
figure(1);
s = tsurf([1 1 1],V,'EdgeColor','none');
c = parula(10);
colormap(c);
colormap summer
% colorbar
shading(gca,'interp')
caxis([0 1])
axis equal
set(gca,'Visible','off');
hold on;

plane_p = [0 X_list{2}(1,2) 0];
plane_n = [0 1 0];
% part behind the plane
[BV,BF] = half_space_intersect(V,F,plane_p,plane_n,'Cap',false);
PI = knnsearch(V,BV);
% but they get glued again...
[BF,I] = cut_edges(BF,sharp_edges(BV,BF));
BV = BV(I,:);
t_back = tsurf(BF,BV,'EdgeColor','none',fsoft);

% part in front of the plane
[FV,FF] = half_space_intersect(V,F,plane_p,-plane_n,'Cap',false);
PI = knnsearch(V,FV);
% but they get glued again...
[FF,I] = cut_edges(FF,sharp_edges(FV,FF));
FV = FV(I,:);
t_front = tsurf(FF,FV,'EdgeColor','none',fsoft);
view(0, 25); 

bg_color = [1.0 1.0 1.0];
camlight;
camproj('persp');
l = light('Position',[4 -4 10],'Style','infinite');
add_shadow(t_back,l,'Color',bg_color*0.8,'BackgroundColor',bg_color,'Fade','infinite');
add_shadow(t_front,l,'Color',bg_color*0.95,'BackgroundColor',bg_color,'Fade','infinite');
apply_ambient_occlusion(t_back,'AddLights',false,'SoftLighting',false,'Factor',0.75);
apply_ambient_occlusion(t_front,'AddLights',false,'SoftLighting',false,'Factor',0.75);
hold off;
for pid = 1:size(parts,1)
  for ii = 1:step
    if ii ~= 2 
      continue; 
    end
    set(s, 'Faces', TF_list{ii}, 'Vertices', X_list{ii}, 'CData', w_list{ii}(:,pid), 'FaceAlpha', 1.0);
    set(s,'facelighting','none');

    % part behind the plane
    plane_p = [0 X_list{ii}(1,2) 0];
    plane_n = [0 1 0];
    % do after cut
    [BV,BF] = half_space_intersect(V,F,plane_p,plane_n,'Cap',false);
    PI = knnsearch(V,BV);
    % but they get glued again...
    [BF,I] = cut_edges(BF,sharp_edges(BV,BF));
    BV = BV(I,:);
    set(t_back, 'Faces',BF, 'Vertices', BV, ...
        'FaceColor',[0 128/255 102/255]*0.6,'FaceAlpha',1.0,'EdgeColor','none',...
         'AmbientStrength', 0.8, 'SpecularStrength', 0.2, 'DiffuseStrength', 0.4);
    
    % part in front of the plane
    [FV,FF] = half_space_intersect(V,F,plane_p,-plane_n,'Cap',false);
    PI = knnsearch(V,FV);
    % but they get glued again...
    [FF,I] = cut_edges(FF,sharp_edges(FV,FF));
    FV = FV(I,:);
    set(t_front, 'Faces', FF, 'Vertices', FV, 'FaceColor',[1.0 1.0 1.0],'FaceAlpha',0.04,'EdgeColor','none',...
       'DiffuseStrength', 0.8);
    
    drawnow;
    % save the png
    print(1, [save_path iges_file(1:end-5) '_part_' num2str(pid) '_step_' num2str(ii) '_cutoff_' num2str(alpha) '.png'], '-dpng', '-r400');
  end
end

%% triangulate the nurbs parts
function [verts,faces,normals] = triangulate_iges(parts)
  faces=[];
  verts=[];
  normals=[];
  for ii=1:numel(parts)
      if isfield(parts{ii}, 'hires_T')
          F = parts{ii}.hires_T;
          V = parts{ii}.hires_x0';
          N = nurbs_normals(parts{ii}.srf.nurbs, parts{ii}.hires_UV, parts{ii}.p);
      else
          F = parts{ii}.T;
          V = parts{ii}.x0';
          N = nurbs_normals(parts{ii}.srf.nurbs, parts{ii}.UV, parts{ii}.p);
      end
      F = F + size(verts,1);
      faces=[faces; F];
      verts=[verts; V];     
      normals=[normals; N];
  end
  % check and remove degenerate faces
  dblA = doublearea(verts,faces);
  deg_fid = find(dblA <= 0);
  faces(deg_fid, :) = [];
end