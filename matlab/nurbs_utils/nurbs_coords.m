function [q, J] = nurbs_coords(srf, u, v)
    % generalized coordinates
    q = zeros(3*size(srf.coefs,2)*size(srf.coefs,3), 1);
    %max_val = max(srf.coefs(1:3,:,:),[],'all');
    for i=1:size(srf.coefs,2)
        for j = 1:size(srf.coefs,3)
            idx = 3*(size(srf.coefs,3)*(i-1) + (j-1)) + 1;
            q(idx:idx+2) = srf.coefs(1:3,i,j);% ./ max_val;
        end
    end
    
    J = zeros(numel(u), numel(v), 3, 3 * prod(srf.number));
    degree = srf.order-1; 
    % Building J matrices (one for each u,v pair)
    for i=1:numel(u)
        si = findspan(srf.number(1)-1, degree(1), u(i), srf.knots{1});
        Ni = basisfun(si, u(i), degree(1), srf.knots{1});
        
        for j=1:numel(v)
            sj = findspan(srf.number(2)-1, degree(2), v(j), srf.knots{2});
            Nj = basisfun(sj, v(j), degree(2), srf.knots{2});
            % outerproduct of B spline values in u & v directions
            % row-wise ordering of matrices in J with u being the row index
            Nij = Ni'*Nj;
            w = srf.coefs(4, si+1-degree(1):si+1, sj+1-degree(2):sj+1);
            
            % Nij is (degree(1)+1 x degree(2)+1) matrix.
            % Bspline decays to zero outside its local support so Nij
            % contains only these nonzero values.
            Nij = Nij .* squeeze(w);
            Nij = Nij ./ sum(Nij,'all');
            for ki=1:size(Nij,1)
                for kj=1:size(Nij,2)
                	Jij = diag(repelem(Nij(ki,kj),3));
                    k = 3*((si-degree(1)+ki-1) * srf.number(2) + (sj-degree(2)+kj-1))+1;
                    J(i,j,:,k:k+2) = Jij;
                end
            end
        end
    end
end
