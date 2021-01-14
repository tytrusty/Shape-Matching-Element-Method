function ME = vem_error_matrix(Y, W, W_S, L, w, E)
    m = size(Y,1);
    ME = zeros(size(L,2), size(L,2));
    d = size(Y,2);
    
    n = size(w,2);
    dk = size(Y,3);

    for i=1:m
%     	Yi = squeeze(Y(i,:,:))*W{i} * W_S{i}; % weighed monomial basis
%         for j=1:n
%             cols=n*dk + d*(j-1)+1 : n*dk + d*j;
%             Yi(:,cols) = eye(d)*w(i,j);
%         end
%         
        I = zeros(d,size(L,2));
        I(:, d*(i-1)+1: d*i) = eye(d);
%         
%         J = I - Yi*L;
        J = I - Y{i}*L;
        
        % Do the multiplication in the loop without L
        % (I-Yi)T(I-Yi) = I - two inner terms + outer term (can add these
        % terms separately, then multiply the L later)
        ME = ME + J'*J;
    end     
end