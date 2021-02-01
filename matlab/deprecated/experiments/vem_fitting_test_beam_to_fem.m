function vem_fitting_test_beam_to_fem
parts = nurbs_from_iges('02_beam2_unbent.igs');
parts_2 = nurbs_from_iges('02_beam2_bent.igs');
figure(1)
clf;
parts=nurbs_plot(parts);

% Finding part correspondences between models using color information. Yea...
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

% Hackiest way in the world to generate centers of mass :)
xcoms = xcoms0;
xcoms(3,:) = xcoms0(3,:) - (0.1)*(xcoms0(1,:).^2.2);
xcoms(1,xcoms(1,:) < 3 & xcoms(1,:) > 2) = 2.2;
xcoms(1,xcoms(1,:) < 4 & xcoms(1,:) > 3) = 2.8;
xcoms(1,xcoms(1,:) > 4) = 3.5;
xcoms(1,xcoms(1,:) < 1) = 0.5;
plot3(xcoms(1,:),xcoms0(2,:),xcoms(3,:),'*','Color','r','MarkerSize',20);

[~, ~, ~, E, x0] = nurbs_assemble_coords(parts);
[~, ~, ~, ~, x2] = nurbs_assemble_coords(parts_2);

[V, ~] = raycast_quadrature(parts, [4 4], 15);
Vplot=plot3(V(1,:),V(2,:),V(3,:),'.','Color','m','MarkerSize',10);
% plot3(x2(1,:),x2(2,:),x2(3,:),'.','Color','r','MarkerSize',10);

order=2; k=9; d=3;

% Compute blending weights
distance_cutoff = 1;

% Blending weights
w = nurbs_blending_weights(parts, V', distance_cutoff);

% Blended undeformed center of mass for old version.
blended_com0 = zeros(size(V,2),3,1);
for i = 1:size(V,2)
    x_com_i = zeros(3,1);
    for j=1:numel(parts)
        x_com_i = x_com_i + w(i,j)*xcoms0(:,j);
    end
    blended_com0(i,:,:) = x_com_i;
end
V0=V;

% L = compute_shape_matrices_newcom(x0, xcoms0, E, [],order, 'hierarchical');
L = compute_shape_matrices(x0, xcoms0, E, order, 'hierarchical');

% x = reshape(J*q,3,[]);

b = [];
for i=1:numel(E)
    b = [b x2(:,E{i}) - xcoms(:,i)];
end
b = b(:);
c = L * b;

for i = 1:size(V,2)
    %%%% NEW VERSION %%%%
    % Dave's suggestion. Separate monomial basis and blending the result
    % of each polynomial.
	V(:,i) = zeros(3,1);
    for j=1:numel(parts)
        Yij = squeeze(monomial_basis_matrix(V0(:,i), xcoms0(:,j), order, k));
        cj = c(d*k*(j-1)+1:d*k*j);
        V(:,i) = V(:,i) + w(i,j) * (Yij * cj + xcoms(:,j));
    end
    %%%%%%%%%%%%%%%%%%%%%

    %%%% OLD VERSION %%%%
    % Created single monomial basis for this point using a blended
    % undeformed center of mass.
    Yi = squeeze(monomial_basis_matrix(V0(:,i), blended_com0(i,:)', order, k));
    blended_com = zeros(3,1); % blended deformed center of mass
    blended_c = zeros(d*k,1); % blended polynomial coefficients
    for j=1:numel(parts)
        blended_com = blended_com + w(i,j)*xcoms(:,j);
        blended_c = blended_c + w(i,j)*c(d*k*(j-1)+1 : d*k*j);
    end
    V(:,i) = Yi*blended_c + blended_com;
    %%%%%%%%%%%%%%%%%%%%%
end

% Plotting
Vplot.XData = V(1,:);
Vplot.YData = V(2,:);
Vplot.ZData = V(3,:);
title('Old approach');
drawnow
end

