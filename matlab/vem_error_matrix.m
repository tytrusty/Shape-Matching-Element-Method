function ME = vem_error_matrix(Y, W, W_S, L, w, E)
    m = size(Y,1);
    ME = zeros(size(L,2), size(L,2));
    d = size(Y,2);
    
    n = size(w,2);
    dk = size(Y,3);
    
%     for ii=1:numel(E)
%        for i=1:numel(E{ii})
%             idx = E{ii}(i);
%             Yi = squeeze(Y(idx,:,:))*W{idx}*W_S{idx};
%             cols=n*dk + d*(ii-1)+1 : n*dk + d*ii;
%             Yi(:,cols) = eye(d);
% 
%             I = zeros(d,size(L,2));
%             I(:, d*(idx-1)+1: d*idx) = eye(d);
% 
%             J = I - Yi*L;
%             ME = ME + J'*J;
%        end
%     end
    
    for i=1:m
    	Yi = squeeze(Y(i,:,:))*W{i} * W_S{i}; % weighed monomial basis
        for j=1:n
            cols=n*dk + d*(j-1)+1 : n*dk + d*j;
            Yi(:,cols) = eye(d)*w(i,j);
        end
        
        I = zeros(d,size(L,2));
        I(:, d*(i-1)+1: d*i) = eye(d);
        
        J = I - Yi*L;
        
        ME = ME + J'*J;
    end     
end