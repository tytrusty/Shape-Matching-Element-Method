function test_cutoff_heuristic

% parts = nurbs_from_iges('puft_no_collar.iges');
% parts = nurbs_from_iges('castle_simple.iges');
parts = nurbs_from_iges('beam2.igs');
figure(1)
clf;
parts=nurbs_plot(parts);

% Assembles global generalized coordinates
[J, ~, q, E, x0] = nurbs_assemble_coords(parts);
    
min_bnd = min(x0, [], 2);
max_bnd = max(x0, [], 2);
max_d = max(max_bnd-min_bnd);

n = 20;
cutoff_vals = linspace(0,max_d, n+1);
cutoff_vals = cutoff_vals(2:end);


enable_secondary_rays = true;
% samples = [5 9 15];
samples = [35 1 1];
[X, ~] = raycast_quadrature(parts, samples(2:3), samples(1));

ratios =zeros(n,1);
nzeros = zeros(n,1);
for i = 1:n
    cutoff = cutoff_vals(i);
    w = distance_weights(parts, X', cutoff, enable_secondary_rays);
    ratio = sum(w(:) > 1e-3) / numel(w)
    w_sum = sum(any(w,2)) / size(w,1)
    ratios(i) = ratio;
    nzeros(i) = w_sum;
end

winner = find(nzeros==1,1,'first');

% plot3(X(1,:),X(2,:),X(3,:),'.','Color','m','MarkerSize',20);

[w, w_I] = nurbs_blending_weights(parts, X', cutoff_vals(winner));
    
[x0_coms, com_cluster, com_map] = generate_com(x0, E, w, numel(parts));
com_plt = plot3(x0_coms(1,:),x0_coms(2,:),x0_coms(3,:), ...
                '.','Color','g','MarkerSize',20);
hold on;
end0