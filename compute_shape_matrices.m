function [B,Q] = compute_shape_matrices(x0, x0_com, E, order)
    
    d = size(x0,1); % dimension

    % Compute number of monomials
    k=0;
    for i=1:order
        k=k+nchoosek(d+i-1,i);
    end

    % TODO only accept cell-types
    if iscell(E)
       Q = cell(numel(E),1);
       B = cell(numel(E),1);
    else
       Q = zeros(k,size(E,2),size(E,1));
       B = zeros(size(E,2),k,size(E,1));
    end
    
    for i=1:size(E,1)
        % Monomial basis of DOF for the i-th shape.
        
        if iscell(E)
            x = x0(:,E{i});
        else
            x = x0(:,E(i,:));
        end
        Qi = monomial_basis(x, x0_com, order);
        
        % Build 'B' matrices (Aqq in shape matching)
        [U, S, V] = svd(Qi*Qi');
        S = diag(S);
        S = 1./ max(S, 1e-4);
        Bi = Qi' * (V * diag(S) * U');

        if iscell(E)
            Q{i} = Qi;
            B{i} = Bi;
        else
        	Q(:,:,i) = Qi;
            B(:,:,i)=Bi;
        end
    end
end

