function Q = monomial_basis(x, x_com, order)
    d = size(x,1);
    
    Q = x - x_com;

    % TODO -- make this general
    if order == 2
        if d==2
            Q(3,:) = Q(1,:).^2;
            Q(4,:) = Q(2,:).^2;
            Q(5,:) = Q(1,:).*Q(2,:);
        else
            Q(5,:) = Q(2,:).^2;
            Q(6,:) = Q(3,:).^2;
            Q(7,:) = Q(4,:).^2;
            Q(8,:) = Q(2,:).*Q(3,:);
            Q(9,:) = Q(3,:).*Q(4,:);
            Q(10,:) = Q(4,:).*Q(2,:);
        end
    end
end

