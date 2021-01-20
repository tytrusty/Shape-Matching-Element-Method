function [x_coms, com_cluster, com_map] = generate_com(parts, x0, E, distance_cutoff, n)
% Generates centers of mass (most of this is temporary)
    centroids = zeros(n, 3);
    x_coms = zeros(3,n);
    for i=1:n
        centroids(i,:) = mean(parts{i}.x0,2)';
    end
    
    cutoff_sqr = distance_cutoff^2;
    adj = zeros(n,n);
    for i=1:n
    	F = parts{i}.hires_T;
        V = parts{i}.hires_x0';
        sqrD = point_mesh_squared_distance(centroids,V,F);
        adj(:,i) = sqrD < cutoff_sqr;
    end
    
    adjacent = cell(n,1);
    for i=1:n
       	adj_idx = find(adj(i,:));
        x_idx = horzcat(E{adj_idx});
        adjacent{i} = x_idx(:);
        x_coms(:,i) = mean(x0(:,adjacent{i}),2);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%
    [C,IA,com_map] = uniquetol(x_coms', 'ByRows', true);
    com_cluster = adjacent(IA);
    x_coms = x_coms(:,IA);
    
end