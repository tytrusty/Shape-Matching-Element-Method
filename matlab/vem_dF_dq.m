function dF_dq = vem_dF_dq(B, dM_dX, E, N, w)
    % It's correct, trust me

    m = size(dM_dX,1);  % # of sampled points
    d = size(dM_dX,3);  % dimension (2 or 3)
    
    dF_dq = zeros(m,d*d,d*N);

    % Rearrange dM_dX to matrix of dM_dX entries stacked column-wise
    dM_dX=permute(dM_dX,[2 3 1]);
    dM_dX=dM_dX(:,:);
    for i = 1:numel(E)
        n = size(B{i},1);

        % gradient contributions from the set E_i
        grad_i = B{i} * dM_dX;
        grad_i = reshape(grad_i,n, d,[]);
        grad_i = bsxfun(@times, grad_i,reshape(w(:,i), [1 1 m]));
        grad_i = permute(grad_i, [3 2 1]);
        idxs = d*(E{i}-1);
        for j = 1:d
            for k = 1:d
                idx=(j-1)*d + k;
                dF_dq(:,idx,idxs+k) = dF_dq(:,idx,idxs+k) + grad_i(:,j,:);
            end
        end
    end

end
