function identify_seams
    % Plot NURBs mesh
    function nurbs = plot_srf(nurbs)
            cm=jet(numel(nurbs));
            for ii = 1:numel(nurbs)
                xi = reshape(nurbs{ii}.x0, 3, nurbs{ii}.subd(1), nurbs{ii}.subd(2));
                plt = surf(squeeze(xi(1,:,:)),squeeze(xi(2,:,:)),squeeze(xi(3,:,:)),'FaceAlpha',0.9,'FaceColor',cm(ii,:));
                hold on;
                nurbs{ii}.plt=plt;
            end
    end

%     part=nurbs_from_iges('rounded_cube.iges',5,0);
    res=repelem(3,16);
    res(5)=6; res(6)=12;
    part=nurbs_from_iges('rocket.iges',res,0);
    part=plot_srf(part);
    axis equal
    
    x0=zeros(3,0);
    E=cell(numel(part),1);
    for i=1:numel(part)
        idx1=size(x0,2)+1;
        x0 = [x0 part{i}.x0];
        idx2=size(x0,2);
        
        % indices into global configuration vector
        part{i}.idx1=idx1;
        part{i}.idx2=idx2;
        E{i}=idx1:idx2;
    end
    
    nurbs=[];
    target=[];
    function f = objective_func(x)
        p = nrbeval(nurbs,x);
        %f = norm(target-p);
        diff = target-p;
        f = diff'*diff;
    end
    
    alpha=1;
    beta=300;
    options = optimoptions('fmincon');
    options.FiniteDifferenceType = 'central';
    options.Display = 'none';

    for i=1:numel(E)
        % Get shape E's u,v values corresponding to each sample.
        [U,V] = meshgrid(part{i}.u, part{i}.v);
        U=U';
        V=V';
        uv = [U(:) V(:)]';
        
        x0_edge = part{i}.x0';
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
        
        nurbs=part{i}.nurbs;
        for j=1:numel(candidates)
            target=x0(:,candidates(j));
            uv_0 = uv0(:,j);
            LB = [part{i}.u_range(1) part{i}.v_range(1)]';
            UB = [part{i}.u_range(2) part{i}.v_range(2)]';
            uv_clamp = min(max(uv_0, LB+[1e-12 1e-12]'), UB-[0.9999999 0.9999999]');
            
            [uv,fval] = fmincon(@objective_func, uv_clamp,[],[],[],[],LB,UB,[],options);
            min_dist(candidates(j)) = min(min_dist(candidates(j)), fval);
        end
        
        dist_threshold = min_dist < alpha;
        dist_threshold(E{i}) = 0;
        verts_to_weld = find(dist_threshold);
        markerColors = hsv(numel(verts_to_weld));
        h1=scatter3(part{i}.x0(1,I(verts_to_weld)),part{i}.x0(2,I(verts_to_weld)),part{i}.x0(3,I(verts_to_weld)),600,'.');
        hold on;
        h2=scatter3(x0(1,verts_to_weld),x0(2,verts_to_weld),x0(3,verts_to_weld),300,'o');
        S_I{i} = verts_to_weld;
        delete(h1);
        delete(h2);
    end
    
end
