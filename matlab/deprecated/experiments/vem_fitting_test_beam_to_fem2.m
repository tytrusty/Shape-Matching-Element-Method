function vem_fitting_test_beam_to_fem2
%polynomial fitting test
%define quadratic projection for 3 point line
%idea: shape match should minimize deformation
%conjecture: corresponsds to keeping higher order terms as close to zero
%as possible
%requirement: surface patch must at least be able to uniquely determine 
% the linear term (deformation gradient)
% this requires explcitly setting the constant term as we've been doing

iges_file = 'beam2.igs';
iges_file = '01_beam_unbent.igs';

% [V,T,F]=readMESH('vem_input.mesh');
% writeOBJ('vem_input.obj',V,F);
% [V,T,F]=readMESH('vem_output.mesh');
% writeOBJ('vem_output.obj',V,F);
% patch('vertices',V,'faces',F,'FaceAlpha',0.2);
% parts_1 = nurbs_from_iges('01_beam_unbent.igs');
% parts_2 = nurbs_from_iges('01_beam_bent.igs');
parts = nurbs_from_iges('02_beam2_unbent.igs');
parts_2 = nurbs_from_iges('02_beam2_bent.igs');
figure(1)
clf;
parts=nurbs_plot(parts);

part_map = zeros(numel(parts),1);
for i=1:numel(parts)
   for j=1:numel(parts_2)
      if parts{i}.srf.color(1) == parts_2{j}.srf.color(1)
          part_map(i) = j;
      end
   end
end
parts_2=parts_2(part_map);
parts_2=nurbs_plot(parts_2);

xcoms = zeros(3, numel(parts));
for i=1:numel(parts)
   xcoms(:,i) = [mean(parts{i}.x0(1,:),2) 0.5 0.5]'; 
end
xcoms(:,1)=xcoms(:,3);
xcoms(:,2)=xcoms(:,end);
xcoms0=xcoms;

xcoms = xcoms0;

[~, ~, q, E, x] = nurbs_assemble_coords(parts);
[~, ~, ~, ~, x2] = nurbs_assemble_coords(parts_2);

% Undeformed Center of mass
x0_com = mean(x,2);

[V, ~] = raycast_quadrature(parts, [4 4], 15);
Vplot=plot3(V(1,:),V(2,:),V(3,:),'.','Color','m','MarkerSize',10);
% plot3(x2(1,:),x2(2,:),x2(3,:),'.','Color','r','MarkerSize',10);

order=2; k=9; d=3;

% Compute blending weights
distance_cutoff = 0.5;
% distance_cutoff = 2;
w = nurbs_blending_weights(parts, V', distance_cutoff);
[W, ~, W_S] = build_weight_matrix(w, d, k, 'Truncate', true);

L = compute_shape_matrices(x, x0_com, E, order, 'hierarchical');
Y = monomial_basis_matrix(V,  x0_com, order, k);

q0 =q;

% x = reshape(J*q,3,[]);

b = [];
for i=1:numel(E)
%     b = [b x2(:,E{i}) - x0_com];
%     b = [b x(:,E{i}) - mean(x2,2)];
    b = [b x2(:,E{i}) - xcoms(:,i)];
end
b = b(:);

c = L * b;
p = c(end-3+1:end);
x_com = x0_com + p; 

for i = 1:size(V,2)
    V(:,i) = squeeze(Y(i,:,:)) * W{i} * W_S{i} * c + x_com_i;
%     V(:,i) = squeeze(Y(i,:,:)) * W{i} * W_S{i} * c + x_com;
end
plot3(x_com(1),x_com(2),x_com(3),'*','Color','r','MarkerSize',20);

% Plotting
Vplot.XData = V(1,:);
Vplot.YData = V(2,:);
Vplot.ZData = V(3,:);
drawnow
end

