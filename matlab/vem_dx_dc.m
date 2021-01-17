function [Y,S] = vem_dx_dc(x, x_coms, w, W_I, com_map, order, k)
    d = size(x,1);
    n = size(x,2);
    
    Y = cell(n,1);
    S = cell(n,1);
    
    for ii=1:n
        m = numel(W_I{ii});
        
        [C,~,ic] = unique(com_map(W_I{ii}));
        ncoms = numel(C);
        
        Yi = zeros(d, d*(k*m + ncoms));
        S_i = zeros(d*(k*m+ncoms), d*(k*size(w,2) + size(x_coms,2)));

        I1 = d*k*m;
        I2 = d*k*size(w,2);
        for i=1:m
            idx = W_I{ii}(i);
            Yij = monomial_basis(x(:,ii), x_coms(:,com_map(idx)), order);
            for j = 1:d
                col_b = d*k*(i-1)+k*(j-1)+1;
                col_e = d*k*(i-1)+k*j;

                Yi(j,col_b:col_e) = w(ii,idx)*Yij;
            end
            I_range = I1+d*(ic(i)-1)+1 : I1+d*ic(i);
            Yi(:,I_range) = Yi(:,I_range) + w(ii,idx)*eye(d);

            % Setting selection matrix terms.
            rows = d*k*(i-1)+1:d*k*i;
            cols = d*k*(idx-1)+1:d*k*idx;
            S_i(rows,cols) = eye(d*k);

            cols = I2+d*(com_map(idx)-1)+1 : I2+d*com_map(idx);
            S_i(I_range, cols) = eye(d);
        end
        Y{ii} = Yi;
        S{ii} = sparse(S_i);
    end
end

