function [w,w_I] = blending_weights_2d(X, V, F, alpha)
    dE = (V(F(:,2),:) - V(F(:,1),:));
    norm2_dE = sum(dE.*dE, 2);

    dVx = X(:,1)' - V((F(:,1)),1);
    dVy = X(:,2)' - V((F(:,1)),2);

    dEdotdV = (dE(:,1).*dVx +  dE(:,2).*dVy)./norm2_dE;

    %bound constraints
    dEdotdV = min(max(dEdotdV, 0),1);

    m=size(X,1);
    w_I = cell(m,1);
    
    w = zeros(size(X,1), size(F,1));

    for jj=1:size(F,1)

        %nearest point on i^th edge of j^th point (column)
        VE_X = (dEdotdV(jj,:).*dE(jj,1) + V(F(jj,1),1))';
        VE_Y = (dEdotdV(jj,:).*dE(jj,2) + V(F(jj,1),2))';

        %build ray from nearest point to i^th point 
        R_X = X(:,1) - VE_X;
        R_Y = X(:,2) - VE_Y;

        dE_X = dE(:,1)';
        dE_Y = dE(:,2)';

        %intersect every edge with very edge 
        det_A = 1./(-R_Y.*dE_X + R_X.*dE_Y);

        %coefficient along the plane
        B_X = X(:,1) - V(F(:,1),1)';
        B_Y = X(:,2) - V(F(:,1),2)';

        first_row = -R_Y.*B_X + R_X.*B_Y;
        second_row = dE_X.*B_Y - dE_Y.*B_X;

        %intersection parameter
        beta = first_row.*det_A;
        gamma = second_row.*det_A;

        beta(gamma < 0) = inf;
        gamma(gamma < 0) = inf;

        gamma(beta > 1) = inf;
        gamma(beta < 0) = inf;

        beta(beta > 1) = inf;
        beta(beta < 0) = inf;

        %deal with special cases
        beta(isinf(det_A)) = inf;
        beta(and(abs(R_X)==0, abs(R_Y) == 0)) = 0;

        gamma(isinf(det_A)) = inf;
        gamma(and(abs(R_X) == 0, abs(R_Y) == 0)) = 0;

        %reconstruct the points I hit
        [~, edge_id] = min(gamma, [], 2);

        hit_points = diag(beta(sub2ind(size(beta), (1:size(X,1))', edge_id)))*dE(edge_id,1:2)  + V(F(edge_id,1),1:2);

        dt = sqrt((hit_points(:,1) - VE_X).^2 + (hit_points(:,2) - VE_Y).^2);
        dx = sqrt((X(:,1) - VE_X).^2 + (X(:,2) - VE_Y).^2);
        w(:,jj) = max(1.0 - dx./min(alpha,dt),0);
    end

    for ii = 1:size(w,1)
        W = diag(1./w(ii,:));
        W(isinf(W)) = 1e8;

        w(ii,:) = quadprog(W, zeros(size(F,1),1), [], [],...
                           ones(1,size(F,1)), 1,...
                           zeros(size(F,1),1), ones(size(F,1),1));
                       
        ToKeep = w(ii,:) > 1e-4;

        % Shape indices for weights above the threshold.
        w_I{ii} = find(ToKeep)';

    end
end