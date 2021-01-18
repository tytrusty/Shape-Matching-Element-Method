function raycasting_integrate_volume
iges_file = 'trident.iges';
% iges_file = 'beam1.igs';

figure(1);
clf;
parts = nurbs_from_iges(iges_file);

cm=jet(numel(parts));
for ii = 1:numel(parts)
    plt = patch('Faces',parts{ii}.hires_T,'Vertices',parts{ii}.hires_x0', ...
        'FaceAlpha',0.3,'EdgeColor','none','FaceColor',cm(ii,:));
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

nsamples = 25;
x_vals = 1:nsamples;
errors = zeros(size(x_vals));

true_vol = 11774.6564; % mm^3 for the trident 
% true_vol = 5; % mm^3 for the beam 

for i = 1:nsamples
   [V, vols] = raycast_quadrature(parts, [i i], i);
   e = (sum(vols) - true_vol)^2;
   errors(i) = e;
   if i == 10
       plot3(V(1,:),V(2,:),V(3,:),'.','Color','m','MarkerSize',20);
   end
end

clf;
semilogy(x_vals,errors);
title('Squared error as the number of samples increase in each direction');
ylabel('Squared Error');
xlabel('# Samples');
end