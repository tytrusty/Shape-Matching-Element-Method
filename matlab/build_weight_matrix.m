function W = build_weight_matrix(w,d,k)
    m = size(w,1); % number of points
    n = size(w,2);  % number of shapes
    
    W = zeros(m, d*d, d*(k*n + 1));
    for i = 1:m
%         Wi = zeros(d*k, d*k*n + d);
        for j = 1:n
            W(i, :, d*k*(j-1)+1:d*k*j) = eye(d*k)*w(i,j);
        end        
    end
end