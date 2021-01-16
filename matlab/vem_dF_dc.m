function [dF_dc, S] = vem_dF_dc(x, x_coms, w, W_I, com_map, order, k)
    d = size(x,1);
    n = size(x,2);

    dF_dc = cell(n,1);
    S = cell(n,1);

    for ii=1:n
        m = numel(W_I{ii});

        dF_dc_i = zeros(d*d, d*k*m);
        S_i = zeros(d*k*m, d*(k*size(w,2) + size(x_coms,2)));
        for jj=1:m
            idx = W_I{ii}(jj);
            dYij = squeeze(monomial_basis_grad(x(:,ii), x_coms(:,com_map(idx)), order));
            
            for i = 1:d
                for j = 1:d
                    row = d*(i-1)+j;
                    col_b = d*k*(jj-1)+k*(j-1)+1;
                    col_e = d*k*(jj-1)+k*j;
                    dF_dc_i(row,col_b:col_e) = w(ii,idx)*dYij(:,i)';
                end
            end

            rows = d*k*(jj-1)+1:d*k*jj;
            cols = d*k*(idx-1)+1:d*k*idx;
            S_i(rows,cols) = eye(d*k);
        end
        dF_dc{ii}=dF_dc_i;
        S{ii} = sparse(S_i);
    end
end

