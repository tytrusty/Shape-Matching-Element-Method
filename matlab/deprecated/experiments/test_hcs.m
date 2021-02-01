function test_hcs
parts = nurbs_from_iges('puft_head.iges');
parts = nurbs_from_iges('beam2.igs');
% parts = nurbs_from_iges('rounded_cube.iges');
figure(1)
clf;
parts=nurbs_plot(parts);

n = numel(parts);
% Assembles global generalized coordinates
[~, ~, ~, E, x0] = nurbs_assemble_coords(parts);

com_threshold = 50;
distance_cutoff = 50;
distance_cutoff = 1.5;

% Generating centers of mass. Temporary method!
x0_coms = generate_com(parts, x0, E, com_threshold, n);

plot3(x0_coms(1,:),x0_coms(2,:),x0_coms(3,:), ...
                '.','Color','r','MarkerSize',20);
hold on;

[V, ~] = raycast_quadrature(parts, [3 3], 10);
Vplot=plot3(V(1,:),V(2,:),V(3,:),'.','Color','m','MarkerSize',10);

[w, w_I] = nurbs_blending_weights(parts, V', distance_cutoff);
% [w, w_I] = nurbs_blending_weights(parts, x0', distance_cutoff);

g = zeros(n,n);  % creat a graph

[wids, ~, ic] = unique(w > 1e-4, 'rows');
w_vals = zeros(size(wids));
w_num = sum(wids,2)
for i=1:size(w,2)
    w_vals(:,i) = accumarray(ic, w(:,i), [], @mean);
end

for i=1:size(wids,1)
    ids = find(wids(i,:));
    if numel(ids) > 1
    	edges = nchoosek(ids,2);
        w_val = w_vals(i,edges(:,1)) + w_vals(i,edges(:,2));
        idx1= sub2ind(size(g),edges(:,1),edges(:,2));
        idx2= sub2ind(size(g),edges(:,2),edges(:,1));
        g(idx1) = g(idx1) + w_val';
        g(idx2) = g(idx2) + w_val';
%         g(idx1) = 1;
%         g(idx2) = 1;
    end
end

G = graph(g);
[mv, mw] = mincut(g)
% H = plot(G,'Layout','layered');
% [mf,asdf,cs,ct] = maxflow(G,n+1,n+2)
% highlight(H,cs,'NodeColor','red')
% highlight(H,ct,'NodeColor','green')
end