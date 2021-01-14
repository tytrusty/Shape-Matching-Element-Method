function dY = monomial_basis_grad_matrix(x, x_coms, w, W_I, com_map, order, k)
    d = size(x,1);
    n = size(x,2);
    
    dY = cell(n,1);
    
    for ii=1:n
        m = numel(W_I{ii});
        
        dYi = zeros(d*d, d*(k*m + size(x_coms,2)));
        for jj=1:m
            idx = W_I{ii}(jj);
            dYij = squeeze(monomial_basis_grad(x(:,ii), x_coms(:,com_map(idx)), order));
            
            for i = 1:d
                for j = 1:d
                    row = d*(i-1)+j;
                    col_b = d*k*(jj-1)+k*(j-1)+1;
                    col_e = d*k*(jj-1)+k*j;
                    dYi(row,col_b:col_e) = w(ii,idx)*dYij(:,i)';
                end
            end
%             
%             for j = 1:d
%                 col_b = d*k*(jj-1)+k*(j-1)+1;
%                 col_e = d*k*(jj-1)+k*j;
%                 Yi(j,col_b:col_e) = w(ii,idx)*Yij;
%             end
        end
        dY{ii}=dYi;
    end
end

