function Y = monomial_basis_matrix(x, x_coms, w, W_I, order, k)
    d = size(x,1);
    n = size(x,2);
    
    Y = cell(n,1);
    
    for ii=1:n
        m = numel(W_I{ii});
        
        Yi = zeros(d, d*(k+1)*m);
        
        I_offset = d*k*m;
        for i=1:m
            idx = W_I{ii}(i);
            Yij = monomial_basis(x(:,ii), x_coms(:,idx), order);
            for j = 1:d
                col_b = d*k*(i-1)+k*(j-1)+1;
                col_e = d*k*(i-1)+k*j;
                Yi(j,col_b:col_e) = w(ii,idx)*Yij;
            end
            I_range = I_offset+d*(i-1)+1 : I_offset+d*i;
            Yi(:,I_range) = w(ii,idx)*eye(d);
        end
        Y{ii}=Yi;
    end
end

