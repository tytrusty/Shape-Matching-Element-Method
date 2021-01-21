function raycasting_integrate_volume
iges_file = 'trident.iges';
% iges_file = 'puft_head.iges';
% iges_file = 'sphere.iges';

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

nsamples = 20;
x_vals = 1:nsamples;
errors = zeros(size(x_vals));

true_vol = 11774.6564; % mm^3 for the trident 
% true_vol = 1.17697114 ; % rounded cube
% true_vol = 5; % mm^3 for the beam 
% true_vol = 267.11
for i = 1:nsamples
   [V, vols] = raycast_quadrature(parts, [i i], 5);
   e = abs(sum(vols) - true_vol) / true_vol;
   errors(i) = e;
   if i == 1
       V_plt = plot3(V(1,:),V(2,:),V(3,:),'.','Color','m','MarkerSize',20);
   else
      V_plt.XData=V(1,:);
      V_plt.YData=V(2,:);
      V_plt.ZData=V(3,:);
      
   end
end

clf;
% semilogy(x_vals,errors);
plot(x_vals,errors);
title('Error as the number of samples increase in each direction');
ylabel('Relative Absolute Error');
xlabel('# Samples');
end