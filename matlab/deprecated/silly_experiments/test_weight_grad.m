function test_weight_grad
%load a mesh
[V,T] = readOBJ([data_dir() '/meshes_obj/star.obj']);
V(1,:) = [];
T = T - 1;
  
F = boundary_faces(T);

function [a,w] = weights(X)
    % From Dave's weight code
    %edge vectors 
    dE = (V(F(:,2),:) - V(F(:,1),:));
    norm2_dE = sum(dE.*dE, 2);

    dVx = X(:,1)' - V((F(:,1)),1);
    dVy = X(:,2)' - V((F(:,1)),2);
    dVz = X(:,3)' - V((F(:,1)),3);

    dEdotdV = (dE(:,1).*dVx +  dE(:,2).*dVy +dE(:,3).*dVz)./norm2_dE;

    %bound constraints
    dEdotdV = min(max(dEdotdV, 0),1);

    a = zeros(size(V,1), size(F,1));
    w = zeros(size(V,1), size(F,1));
    
    %%% THIS BIT DOES THE RAYTRACING WEIGHTS %%%%%
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

        %beta(:,jj) = inf(size(V,1),1);
        %gamma(:,jj) = inf(size(V,1),1);

        gamma(beta > 1) = inf;
        gamma(beta < 0) = inf;

        beta(beta > 1) = inf;
        beta(beta < 0) = inf;

        %deal with special cases
        beta(isinf(det_A)) = inf;
        beta(and(abs(R_X) ==0, abs(R_Y) == 0)) = 0;

        gamma(isinf(det_A)) = inf;
        gamma(and(abs(R_X) == 0, abs(R_Y) == 0)) = 0;

        %reconstruct the points I hit
        [min_gamma, edge_id] = min(gamma, [], 2);

        hit_points = diag(beta(sub2ind(size(beta), (1:size(V,1))', edge_id)))*dE(edge_id,:)  + V(F(edge_id,1),:);

        dt = sqrt((hit_points(:,1) - VE_X).^2 + (hit_points(:,2) - VE_Y).^2);
        dx = sqrt((X(:,1) - VE_X).^2 + (X(:,2) - VE_Y).^2);
        w(:,jj) = max(1.0 - dx./0.5,0);
%         w(:,jj) = max(1.0 - dx./min(0.5,dt),0);
    end

    for ii = 1:size(a,1)
        
        W = diag(1./w(ii,:));
        W(isinf(W)) = 1e8;
        
        a(ii,:) = quadprog(W, zeros(size(F,1),1), [], [],...
                            ones(1,size(F,1)), 1,...
                           zeros(size(F,1),1), ones(size(F,1),1));
        
    end
%     w = a;
end

h = 1e-8;
[w,a]= weights(V);
w_posx = weights(V + [h 0 0]);
w_negx = weights(V - [h 0 0]);
w_posy = weights(V + [0 h 0]);
w_negy = weights(V - [0 h 0]);
dw_dx = (w_posx - w_negx) / (2 * h);
dw_dy = (w_posy - w_negy) / (2 * h);
find(max(dw_dx(:,1)))
find(max(dw_dy(:,1)))
max(dw_dx(:,1))
max(dw_dy(:,1))
C=rescale(dw_dx);

%make a little UI to play with
f = figure(1);
clf;
hold on
% srf = tsurf(T, V);
map = parula(256);
w_i = floor(a(:,1)*255) + 1;
C = map(w_i,:);
p=patch('Faces',T,'Vertices',V, 'FaceColor', 'interp', 'EdgeColor', 'none');
p.FaceVertexCData = C;
p.CDataMapping = 'scaled';
quiver(V(:,1),V(:,2),dw_dx(:,1), dw_dy(:,1), 'r');
axis equal

drawnow;
hold off

lengths = sqrt(dw_dx(:,1).*dw_dx(:,1) + dw_dy(:,1).*dw_dy(:,1));
end
