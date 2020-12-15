function S_I = identify_seams(parts, x0, E, alpha, beta)
    if nargin < 5
        beta=300;
    end
    
    if nargin < 4
       alpha=1; 
    end

    nurbs=[];
    target=[];
    function f = objective(x)
        p = nrbeval(nurbs,x);
        diff = target-p;
        f = diff'*diff;
    end

    options = optimoptions('fmincon');
    options.FiniteDifferenceType = 'central';
    options.Display = 'none';
    
    S_I = cell(numel(E),1);
    
    tic
    for i=1:numel(E)
        % Get shape E's u,v values corresponding to each sample.
        [U,V] = meshgrid(parts{i}.u, parts{i}.v);
        U=U';
        V=V';
        uv = [U(:) V(:)]';
        
        x0_edge = parts{i}.x0';
        diff_x = x0(1,:) - x0_edge(:,1);
        diff_y = x0(2,:) - x0_edge(:,2);
        diff_z = x0(3,:) - x0_edge(:,3);
        diff = diff_x.^2 + diff_y.^2 + diff_z.^2;
        [min_dist,I] = min(diff,[],1);
        dist_threshold = min_dist < beta;
        dist_threshold(E{i}) = 0;
        
        % Candidate points for which we're optimizing
        candidates = find(dist_threshold);
        uv0 = uv(:, I(candidates)); % initial uv values
        
        nurbs=parts{i}.nurbs;
        for j=1:numel(candidates)
            if min_dist(candidates(j)) > alpha
                target=x0(:,candidates(j));
                uv_0 = uv0(:,j);
                LB = [parts{i}.u_range(1) parts{i}.v_range(1)]';
                UB = [parts{i}.u_range(2) parts{i}.v_range(2)]';
                uv_clamp = min(max(uv_0, LB+[1e-12 1e-12]'), UB-[0.9999999 0.9999999]');

                [uv,fval] = fmincon(@objective, uv_clamp,[],[],[],[],LB,UB,[],options);
                min_dist(candidates(j)) = min(min_dist(candidates(j)), fval);
            end
        end
        
        dist_threshold = min_dist < alpha;
        dist_threshold(E{i}) = 0;
        verts_to_weld = find(dist_threshold);
        S_I{i} = verts_to_weld;
    end
    toc
    
end
