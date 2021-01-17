function test_raycast_mug
iges_file = 'mug_models/mug_intersect_2.igs';
% iges_file = 'beams/beam2.igs';
iges_file = 'starship.iges';
% iges_file = 'trident.iges';
% iges_file = 'wrench.iges';

fig=figure(1);
clf;
parts = nurbs_from_iges(iges_file);

cm=jet(numel(parts));
for ii = 1:numel(parts)
    plt = patch('Faces',parts{ii}.hires_T,'Vertices',parts{ii}.hires_x0', ...
        'FaceAlpha',0.2,'EdgeColor','none','FaceColor',cm(ii,:));
    hold on;
    parts{ii}.plt=plt;
end
    
set(gcf,'color','w');
axis equal
lighting gouraud;
lightangle(gca,0, 40)
xlabel('x');
zlabel('z');
view(-15,30);


[V, ~] = raycast_quadrature(parts, [1 20], 5);
[V, ~] = raycast_quadrature(parts, [3 3], 5);
plot3(V(1,:),V(2,:),V(3,:),'.','Color','m','MarkerSize',20);
end