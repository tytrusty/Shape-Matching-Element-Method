clear;
iges_file = 'puft_no_collar.iges';
parts = nurbs_from_iges(iges_file);
[V,F,N,C] = triangulate_iges(parts);

save_path = ['../output/'];
mkdir(save_path); 
clf;
figure(1);
t = tsurf(F,V,'CData', C, 'VertexNormals', N, 'EdgeColor', 'none',...
          fsoft, 'FaceAlpha', 0.9, ...
          'AmbientStrength', 0.8, 'SpecularStrength', 0.2, 'DiffuseStrength', 0.4);
hold on
  CM = cbrewer('RdYlBu',20);
  colormap(CM);
shading(gca,'interp')
axis equal
set(gca,'Visible','off');
bg_color = [1.0 1.0 1.0];
set(gcf,'Color',bg_color);
camlight;
camproj('persp');
view(0, 25);
l = light('Position',[4 -4 10],'Style','infinite');
add_shadow(t,l,'Color',bg_color*0.9,'BackgroundColor',bg_color,'Fade','infinite');
apply_ambient_occlusion(t,'AddLights',false,'SoftLighting',false,'Factor',0.75);
hold off

print(1, [save_path iges_file(1:end-5) '.png'], '-dpng', '-r400');

%% triangulate the nurbs parts
function [verts,faces,normals,colors] = triangulate_iges(parts)
  faces=[];
  verts=[];
  normals=[];
  colors=[];
 
  perm = randperm(numel(parts));
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
      colors = [colors; repmat(perm(ii), size(V,1), 1)];
  end
  % check and remove degenerate faces
  dblA = doublearea(verts,faces);
  deg_fid = find(dblA <= 0);
  faces(deg_fid, :) = [];
end
