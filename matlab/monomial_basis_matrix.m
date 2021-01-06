function Y = monomial_basis_matrix(x, x_com, order, k)
    d = size(x,1);
    n = size(x,2);
    
    Y = zeros(n, d, d*k);
    
    Q = x - x_com;
    

    % TODO -- make this general
    if order == 2
        if d==2
            Q(3,:) = Q(1,:).^2;
            Q(4,:) = Q(2,:).^2;
            Q(5,:) = Q(1,:).*Q(2,:);
        else
            Q(4,:) = Q(1,:).^2;
            Q(5,:) = Q(2,:).^2;
            Q(6,:) = Q(3,:).^2;
            Q(7,:) = Q(1,:).*Q(2,:);
            Q(8,:) = Q(2,:).*Q(3,:);
            Q(9,:) = Q(3,:).*Q(1,:);
        end
    end
    Q = Q';
    
    for i=1:d
        col_range = k*(i-1)+1 : k*i;
    	Y(:, i, col_range) = Q;
    end
end

