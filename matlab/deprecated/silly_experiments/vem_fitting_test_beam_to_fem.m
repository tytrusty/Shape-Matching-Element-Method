function vem_fitting_test_beam_to_fem
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

xcoms(3,:) = xcoms0(3,:) - (0.1)*(xcoms0(1,:).^2.2);
xcoms(1,xcoms(1,:) < 3 & xcoms(1,:) > 2) = 2.2;
xcoms(1,xcoms(1,:) < 4 & xcoms(1,:) > 3) = 2.8;
xcoms(1,xcoms(1,:) > 4) = 3.5;
xcoms(1,xcoms(1,:) < 1) = 0.5;
COMplot=plot3(xcoms(1,:),xcoms0(2,:),xcoms(3,:),'*','Color','r','MarkerSize',20);

[~, ~, q, E, x0] = nurbs_assemble_coords(parts);
[~, ~, ~, ~, x2] = nurbs_assemble_coords(parts_2);

% Undeformed Center of mass
x0_com = mean(x0,2);

[V, ~] = raycast_quadrature(parts, [4 4], 15);
Vplot=plot3(V(1,:),V(2,:),V(3,:),'.','Color','m','MarkerSize',10);
% plot3(x2(1,:),x2(2,:),x2(3,:),'.','Color','r','MarkerSize',10);

order=2; k=9; d=3;

% Compute blending weights
distance_cutoff = 0.5;
% distance_cutoff = 2;
w = nurbs_blending_weights(parts, V', distance_cutoff);
[W, ~, W_S] = build_weight_matrix(w, d, k, 'Truncate', true);

% L = compute_shape_matrices(x0, x0_com, E, order, 'hierarchical');
% Y = monomial_basis_matrix(V,  x0_com, order, k);

x0_com_i = zeros(size(V,2),3,1);
for i = 1:size(V,2)
    x_com_i = zeros(3,1);
    for j=1:numel(parts)
        x_com_i = x_com_i + w(i,j)*xcoms0(:,j);
    end
    x0_com_i(i,:,:) = x_com_i;
end


L = compute_shape_matrices(x0, xcoms0, E, order, 'hierarchical');
Y = monomial_basis_matrix(V,  x0_com_i', order, k);

q0 =q;

% x = reshape(J*q,3,[]);

b = [];
for i=1:numel(E)
%     b = [b x(:,E{i}) - x0_com];
%     b = [b x(:,E{i}) - mean(x,2)];
    b = [b x2(:,E{i}) - xcoms(:,i)];
%     b = [b x(:,E{i}) - xcoms0(:,i)];
%     b = [b x(:,E{i}) ];
end
b = b(:);

c = L * b;
p = c(end-3+1:end);
x_com = x0_com + p; 

for i = 1:size(V,2)
    x_com_i = zeros(3,1);
    for j=1:numel(parts)
        x_com_i = x_com_i + w(i,j)*xcoms(:,j);
    end
    V(:,i) = squeeze(Y(i,:,:)) * W{i} * W_S{i} * c + x_com_i;
%     V(:,i) = squeeze(Y(i,:,:)) * W{i} * W_S{i} * c + x_com;
end

% Plotting
Vplot.XData = V(1,:);
Vplot.YData = V(2,:);
Vplot.ZData = V(3,:);
drawnow
end

