function vem_fitting_test_beam
%polynomial fitting test
%define quadratic projection for 3 point line
%idea: shape match should minimize deformation
%conjecture: corresponsds to keeping higher order terms as close to zero
%as possible
%requirement: surface patch must at least be able to uniquely determine 
% the linear term (deformation gradient)
% this requires explcitly setting the constant term as we've been doing

iges_file = 'beam2.igs';
parts = nurbs_from_iges(iges_file);
figure(1)
clf;
parts=nurbs_plot(parts);

xcoms = zeros(3, numel(parts));
for i=1:numel(parts)
   xcoms(:,i) = [mean(parts{i}.x0(1,:),2) 0.5 0.5]'; 
end
xcoms(:,end-1)=xcoms(:,1);
xcoms(:,end)=xcoms(:,end-2);
COMplot=plot3(xcoms(1,:),xcoms(2,:),xcoms(3,:),'*','Color','r','MarkerSize',10);

xcoms0=xcoms;

[J, ~, q, E, x0] = nurbs_assemble_coords(parts);

% Undeformed Center of mass
x0_com = mean(x0,2);

[V, ~] = raycast_quadrature(parts, [3 3], 15);
Vplot=plot3(V(1,:),V(2,:),V(3,:),'.','Color','m','MarkerSize',10);

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
        x_com_i = x_com_i + w(i,j)*xcoms(:,j);
    end
    x0_com_i(i,:,:) = x_com_i;
end


L = compute_shape_matrices(x0, xcoms0, E, order, 'hierarchical');
Y = monomial_basis_matrix(V,  x0_com_i', order, k);

q0 =q;

steps=100;
for jj = 1:steps
alpha=jj/steps;
beta=35;
q(3:3:end) = q0(3:3:end) - (alpha/beta)*(q0(1:3:end).^3);
xcoms(3,:) = xcoms0(3,:) - (alpha/beta)*(xcoms0(1,:).^3);
x = reshape(J*q,3,[]);

b = [];
for i=1:numel(E)
%     b = [b x(:,E{i}) - x0_com];
%     b = [b x(:,E{i}) - mean(x,2)];
    b = [b x(:,E{i}) - xcoms(:,i)];
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
COMplot.ZData = xcoms(3,:);
x_idx=0;
for i=1:numel(parts)
    x_sz = size(parts{i}.x0,2);
    xi = x(:,x_idx+1:x_idx+x_sz);
    parts{i}.plt.Vertices =xi';
    x_idx = x_idx+x_sz;
end
% plot3(q(1:3:end),q(2:3:end),q(3:3:end),'*','Color','g');
pause(0.05)
drawnow

end
end

