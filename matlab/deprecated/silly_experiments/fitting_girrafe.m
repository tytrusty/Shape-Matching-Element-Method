function fitting_girrafe
parts = nurbs_from_iges('sword.igs');
figure(1)
clf;
parts=nurbs_plot(parts);

n = numel(parts);
% Assembles global generalized coordinates
[J, hires_J, q, E, x0] = nurbs_assemble_coords(parts);

com_threshold =80;
% Generating centers of mass. Temporary method!
[x0_coms, com_cluster, com_map] = generate_com(parts, x0, E, ...
    com_threshold, n);

com_plt = plot3(x0_coms(1,:),x0_coms(2,:),x0_coms(3,:), ...
                '.','Color','r','MarkerSize',20);
hold on;

[V, ~] = raycast_quadrature(parts, [8 8], 15);
Vplot=plot3(V(1,:),V(2,:),V(3,:),'.','Color','m','MarkerSize',10);
% plot3(x2(1,:),x2(2,:),x2(3,:),'.','Color','r','MarkerSize',10);

order=2; k=9; d=3;

% Compute blending weights
distance_cutoff = 50;

% x = reshape(J*q,3,[]);
q(3:3:end) = q(3:3:end) + (1/500)*(q(2:3:end).^2);
x2 = reshape(J*q,3,[]);

x = reshape(hires_J*q,3,[]);

% Update NURBs plots
x_idx=0;
for i=1:numel(parts)
    x_sz = size(parts{i}.hires_x0,2);
    xi = x(:,x_idx+1:x_idx+x_sz);
    parts{i}.plt.Vertices =xi';
    x_idx = x_idx+x_sz;
end


% Compute Shape weights
[w, w_I] = nurbs_blending_weights(parts, V', distance_cutoff);
[w0, w0_I] = nurbs_blending_weights(parts, x0', distance_cutoff);

% Build Monomial bases for all quadrature points
[Y,Y_S] = vem_dx_dc(V, x0_coms, w, w_I, com_map, order, k);
[Y0,Y0_S] = vem_dx_dc(x0, x0_coms, w0, w0_I, com_map, order, k);

V0=V;

L = compute_shape_matrices(x0, x0_coms, com_map, E, com_cluster, ...
    order, 'hierarchical');

b = [];
for i=1:numel(E)
    b = [b x2(:,E{i})];
end
b = b(:);
c = L * b;

for i = 1:size(V,2)
	V(:,i) = zeros(3,1);
    for j=1:numel(parts)
        V(:,i) = Y{i} * Y_S{i} * c;
    end
end

% Plotting
Vplot.XData = V(1,:);
Vplot.YData = V(2,:);
Vplot.ZData = V(3,:);
title('Old approach');
drawnow
end

