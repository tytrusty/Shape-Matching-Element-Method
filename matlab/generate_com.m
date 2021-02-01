function [x_coms, com_cluster, com_map] = generate_com(x0, E, w, n)
    one_com = false; % set true if you just want one com %todo make param
    
    if (~one_com)
        % Find all unique center of mass candidates indicated by weights
        % with the same code.
        % "Code" refers to the bitstring corresponding to nonzero weights
        %   e.g. w=[0.1 0.2 0 0 0.7] would be a code of wid=[1 1 0 0 1]
        [wids, ~, ic] = unique(w > 1e-4, 'rows');
        w_vals = zeros(size(wids)); % edge weights for each code -> parts
        w_num = sum(wids,2);        % # of parts per code
        
        % Setting weights of each 
        for i=1:size(w,2)
            w_vals(:,i) = accumarray(ic, w(:,i), [], @sum);
        end
        
        % Remove candidates consisting of a single part
        singletons = (w_num==1);
        w_vals(singletons,:) = [];
        wids(singletons,:) = [];
        w_num(singletons) = [];

        m = size(wids,1);               % number of candidate coms
        ninternal = sum(w_num);         % # of candidate -> part edges
        nedges = ninternal + m + n;     % total # of edges
        
        % Equality constraints that satisfy flow conservation constraint.
        Aeq = zeros(m + n, nedges);     % # nodes x # edges
        beq = zeros(m+n,1);
        
        % Edge map maps com candidate --> part edges to positions in the
        % constraint matrix.
        edge_map = zeros(size(wids));  
        
        % Cost assigned to each edge for min-cost minimization.
        cost_f = ones(nedges,1);
        
        % Add com conservation constraints
        idx=0;
        for i=1:m
            range = idx+1:idx+w_num(i);
            % candidate inflow == sum of candidate outflow
            Aeq(i, ninternal + i) = 1;
            Aeq(i, range) = -1;

            map_loc = wids(i,:);
            edge_map(i, map_loc) = range;
            
            % Negativing w_vals so that high weight sums result in largest
            % minimization of the 
            cost_f(range) = - w_vals(i, map_loc);
            idx=idx+w_num(i);
        end

        % Add part conservation constraints. Enforce that all the sum
        % of the flow to each part is equal to the parts' outgoing flows to
        % the sink, t.
        for i=1:n
            Aeq(m + i, ninternal + m + i) = 1;
            Aeq(m + i, nonzeros(edge_map(:,i))) = -1;
        end
        
        % Setting flow capacity constraints
        lb=zeros(nedges,1);
        lb(ninternal + m + 1:end) = 1;       % Guarantee each part has com
        ub=ones(nedges,1);
        ub(ninternal+1:ninternal+m) = w_num; % Capacity on source edges
        
        % If we were just doing max flow (without costs):
        % f = zeros(nedges,1);
        % f(ninternal+1:ninternal+m)= -1;
        
        % ILP with integer constraints on all flow values.
        x = intlinprog(cost_f, ones(nedges,1),[],[],Aeq,beq,lb,ub);
        
        % Unused, but uncomment to verify correctness.
        % s_flows = x(ninternal+1:ninternal+m);
        % t_flows = x(ninternal+m+1:end);
        
        % Assign centers of mass to edges
        com_map = zeros(n,1);
        for i = 1:m
            [~,j,v] = find(edge_map(i,:));
            com_map(j(x(v) > 0)) =i;
        end
        
        active_coms = unique(com_map);
        com_cluster = cell(numel(active_coms), 1);
        x_coms = zeros(3,numel(active_coms));

        for i = 1:numel(active_coms)
            % Find all parts use in computing this center of mass and
            % form list of all points used to compute the mean.
            part_idx = wids(active_coms(i),:);
            com_cluster{i} = vertcat(E{part_idx});
            
            % Output initial (undeformed) center of mass position.
            x_coms(:,i) = mean(x0(:,com_cluster{i}),2);

            % Find parts that use this center of mass as their element's
            % center of mass and update their com_map index.
            com_map(com_map == active_coms(i)) = i;
        end
    else
        x_coms = mean(x0,2);
        com_map = ones(n,1);
        com_cluster = {1:size(x0,2)};
    end
end