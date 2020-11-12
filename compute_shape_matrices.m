function [B,Q] = compute_shape_matrices(x0, x0_com, E, order)
    if order == 2
        k = 5;
    else
        k = 2;
    end
    
    Q = zeros(k,size(E,2),size(E,1));
    B = zeros(size(E,2),k,size(E,1));
    
    for i=1:size(E,1)
        % Monomial basis of DOF for the i-th shape.
        Qi = x0(:,E(i,:)) - x0_com;
        if order == 2
        	Qi(3,:) = Qi(1,:).^2;
            Qi(4,:) = Qi(2,:).^2;
            Qi(5,:) = Qi(1,:).*Qi(2,:);
        end
        Q(:,:,i) = Qi;
        
        % Build 'B' matrices (Aqq in shape matching)
        [U, S, V] = svd(Qi*Qi');
        S = diag(S);
        S = 1./ max(S, 1e-4);
        Bi = Qi' * (V * diag(S) * U');
        B(:,:,i)=Bi;
    end

end

