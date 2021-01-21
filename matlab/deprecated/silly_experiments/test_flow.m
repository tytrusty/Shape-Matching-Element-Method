function test_flow
% parts = nurbs_from_iges('mug_complex.iges');
parts = nurbs_from_iges('gumby.iges');
% parts = nurbs_from_iges('beam2.igs');
% parts = nurbs_from_iges('T.iges');
% parts = nurbs_from_iges('rounded_cube.iges');
figure(1)
clf;
parts=nurbs_plot(parts);

n = numel(parts);
% Assembles global generalized coordinates
[~, ~, ~, E, x0] = nurbs_assemble_coords(parts);

% distance_cutoff = 50; % puft_simple
% distance_cutoff = 8; % for T
% distance_cutoff = 1; % beam
distance_cutoff = 101; % beam

% Generating centers of mass. Temporary method!
% x0_coms = generate_com(parts, x0, E, com_threshold, n);
% 
% plot3(x0_coms(1,:),x0_coms(2,:),x0_coms(3,:), ...
%                 '.','Color','r','MarkerSize',20);
% hold on;

[V, ~] = raycast_quadrature(parts, [11 8], 25);
Vplot=plot3(V(1,:),V(2,:),V(3,:),'.','Color','m','MarkerSize',10);

[w, w_I] = nurbs_blending_weights(parts, V', distance_cutoff);
% [w, w_I] = nurbs_blending_weights(parts, x0', distance_cutoff);

[wids, ~, ic] = unique(w > 1e-4, 'rows');
w_vals = zeros(size(wids));
w_num = sum(wids,2);
for i=1:size(w,2)
    w_vals(:,i) = accumarray(ic, w(:,i), [], @sum);
end

% Remove single coms
w_vals(w_num==1,:) = [];
wids(w_num==1,:) = [];
w_num(w_num==1) = [];

m=size(wids,1); % number of candidate coms
ninternal = sum(w_num);
nedges = ninternal + m + n; % total # of edges

Aeq = zeros(m + n, nedges); % # nodes x # edges
beq = zeros(m+n,1);
edge_map = zeros(size(wids));

cost_f = ones(nedges,1);

% Add com conservation constraints
idx=0;
for i=1:m
    range = idx+1:idx+w_num(i);
    Aeq(i, ninternal + i) = 1;
    Aeq(i, range) = -1;
    
    map_loc = wids(i,:);
    edge_map(i, map_loc) = range;
    cost_f(range) = - w_vals(i, map_loc);
    idx=idx+w_num(i);
end

% Add part conservation constraints
for i=1:n
    Aeq(m + i, ninternal + m + i) = 1;
    Aeq(m + i, nonzeros(edge_map(:,i))) = -1;
end

lb=zeros(nedges,1);
lb(ninternal + m + 1:end) = 1; % a flow must occur!
ub=ones(nedges,1);
ub(ninternal+1:ninternal+m)= w_num; % source edges have higher capacity

% Just max flow (without costs)
f = zeros(nedges,1);
f(ninternal+1:ninternal+m)= -1;

% x = intlinprog(f, ones(nedges,1),[],[],Aeq,beq,lb,ub);
x = intlinprog(cost_f, ones(nedges,1),[],[],Aeq,beq,lb,ub);

s_vals = x(ninternal+1:ninternal+m);
t_vals = x(ninternal+m+1:end);

com_map = zeros(n,1);
for i = 1:m
    [~,j,v] = find(edge_map(i,:));
    com_map(j(x(v) > 0)) =i;
end

active_coms = unique(com_map);
for i = 1:numel(active_coms)
    nparts = parts(wids(active_coms(i),:));
    x_com = zeros(3,1);
    num = 0;
    for j=1:numel(nparts)
        num = num + size(nparts{j}.x0,2);
        x_com = x_com + sum(nparts{j}.x0,2);
    end
    x_com = x_com ./ num;
    plot3(x_com(1,:),x_com(2,:),x_com(3,:), ...
                '.','Color','r','MarkerSize',20);
    hold on;
    
    assoc_parts = find(com_map == active_coms(i))
    for j=1:numel(assoc_parts)
       centroid = mean(parts{assoc_parts(j)}.x0,2);
       plot3([centroid(1); x_com(1)],[centroid(2); x_com(2)],[centroid(3); x_com(3)], ...
                '-','Color','r','MarkerSize',20);
    end
end