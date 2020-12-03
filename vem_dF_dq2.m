function dF_dq = vem_dF_dq2(B, dM_dX, E, N, w)
    m = size(dM_dX,1);  % # of sampled points
    d = size(dM_dX,3);  % dimension (2 or 3)
    
    dF_dq = zeros(m,d*d,d*N);
    
    % for i = 1:|E|
    %   Bi
    %   n = size(B,1)
    %   dcomdq = -1/N , n
    %   B * [ dM1x | dM1y | dM2x | dM2y | ... ]
    
    % Rearrange dM_dX to matrix of dM_dX entries stacked column-wise
    
    dM_dX=permute(dM_dX,[2 3 1]);
    dM_dX=dM_dX(:,:);
    for i = 1:numel(E)
        n = size(B{i},1);
        
        % gradient of center of mass w.r.t q
        dcom_dq = repelem(-1/N,n);
        grad_com = dcom_dq * B{i} * dM_dX;
        grad_com = reshape(grad_com,d,[])';
        grad_com = bsxfun(@times, grad_com,w(:,i));
        grad_com = repmat(grad_com,[1 1 N]);     
        
        % gradient contributions from the set E_i
        grad_i = (1 - 1/N) * B{i} * dM_dX;
        grad_i = reshape(grad_i,n, d,[]);
        grad_i = bsxfun(@times, grad_i,reshape(w(:,i), [1 1 N]));
        grad_i = permute(grad_i, [3 2 1]);
        idxs = d*(E{i}-1);
        for j = 1:d
            for k = 1:d
                idx=(j-1)*d + k;
                dF_dq(:,idx,k:d:end) = dF_dq(:,idx,k:d:end) + grad_com(:,j,:);
                dF_dq(:,idx,idxs+k) = dF_dq(:,idx,idxs+k) + grad_i(:,j,:);
                size(dF_dq(:,idx,idxs+k))
                %dF_dq(:,idx,) = (N-1)*grad(:,j,:);
            end
        end
%         grad=grad*(N-1);
%         
%         for j = 1:numel(E{i})
%             for k = 1:d
%                 idx = d*(E(j)-1);
%                 J(k, idx + k, i) = squeeze(J(k, idx + k, i)) + BM_i(j);
%             end
%         end
    end

end
