function [f, f_grad, g, g_grad] = nurbs_grad(srf, u, v)

    % Computing gradient w.r.t to each u,v pair.
    f = zeros(numel(u), numel(v), 3, 3 * prod(srf.number));
    f_grad = zeros(numel(u), numel(v), 2, 3, 3 * prod(srf.number));
    g = zeros(numel(u), numel(v), 1);
    g_grad = zeros(numel(u), numel(v), 2);
    
    degree = srf.order-1; 
    % Building J matrices (one for each u,v pair)
    for i=1:numel(u)
        si = findspan(srf.number(1)-1, degree(1), u(i), srf.knots{1});
        Ni = basisfun(si, u(i), degree(1), srf.knots{1});
        % Input: i,pl,ml   %knot span, degree of curve, number of control points
        %        u , u_knotl  %parameter value, vector of knots
        Ni_grad = dersbasisfuns(si,degree(1),u(i),srf.knots{1},1);
        
        for j=1:numel(v)
            sj = findspan(srf.number(2)-1, degree(2), v(j), srf.knots{2});
            Nj = basisfun(sj, v(j), degree(2), srf.knots{2});
            Nj_grad = dersbasisfuns(sj,degree(2),v(j),srf.knots{2},1);
            
            % outerproduct of B spline values in u & v directions
            % row-wise ordering of matrices in J with u being the row index
            Nij = Ni'*Nj;
            Nij_ugrad = Ni_grad'*Nj;
            Nij_vgrad = Ni'*Nj_grad;
            
            w = srf.coefs(4, si+1-degree(1):si+1, sj+1-degree(2):sj+1);
            
            % Nij is (degree(1)+1 x degree(2)+1) matrix.
            % Bspline decays to zero outside its local support so Nij
            % contains only these nonzero values.
            Nij = Nij .* squeeze(w);
            Nij_ugrad = Nij_ugrad .* squeeze(w);
            Nij_vgrad = Nij_vgrad .* squeeze(w);

            % denominator terms
            g(i,j) = sum(Nij,'all'); 
            g_grad(i,j,1) = sum(Nij_ugrad,'all'); 
            g_grad(i,j,2) = sum(Nij_vgrad,'all'); 
            
            for ki=1:size(Nij,1)
                for kj=1:size(Nij,2)
                    k = 3*((si-degree(1)+ki-1) * srf.number(2) + (sj-degree(2)+kj-1))+1;
                    Jij = diag(repelem(Nij(ki,kj),3));
                    Jij_ugrad = diag(repelem(Nij_ugrad(ki,kj),3));
                    Jij_vgrad = diag(repelem(Nij_vgrad(ki,kj),3));
                    f(i,j,:,k:k+2) = Jij;
                    f_grad(i,j,1,:,k:k+2) = Jij_ugrad;
                    f_grad(i,j,2,:,k:k+2) = Jij_vgrad;
                end
            end
        end
    end
end