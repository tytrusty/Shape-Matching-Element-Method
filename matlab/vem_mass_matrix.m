function M = vem_mass_matrix(B, Q, w, d, N, E)
    m = size(Q,2);
    M = zeros(d*(N+1), d*(N+1));
    J = zeros(d,d*N,m);
    for i = 1:size(E,1)
        n=size(B{i},1);
        w_i = reshape(w(:,i), [1 1 m]);
        J = J + bsxfun(@times, vem_jacobian(B{i},Q,n,d,N,E{i}), w_i);
    end
    
    P = eye(d);
    P = repmat(P,N,1);
    for j=1:m        
%         % Top Left
%         M(1:d*N,1:d*N) = M(1:d*N,1:d*N) + J(:,:,j)'*J(:,:,j);
%         
%         % Bottom Left
%         M(d*N + 1:d*(N+1),1:d*N) = M(d*N + 1:d*(N+1),1:d*N) + J(:,:,j);
%         
%         % Top Right
%         M(1:d*N, d*N + 1:d*(N+1)) = M(1:d*N, d*N + 1:d*(N+1)) + J(:,:,j)';
        % Top Left
%         JTJ = J(:,:,j)'*J(:,:,j);
%         
%         M(1:d*N,1:d*N) = M(1:d*N,1:d*N) + JTJ;
%         
%         % Bottom Left
%         M(d*N + 1:d*(N+1),1:d*N) = M(d*N + 1:d*(N+1),1:d*N) - P' * JTJ;
%         
%         % Top Right
%         M(1:d*N, d*N + 1:d*(N+1)) = M(1:d*N, d*N + 1:d*(N+1)) - JTJ * P;
%         
%         % Bottom right block
%         M(d*N+1:d*(N+1), d*N+1:d*(N+1)) = M(d*N+1:d*(N+1), d*N+1:d*(N+1)) + P'*JTJ*P;

        JTJ = J(:,:,j)'*J(:,:,j);
        
        M(1:d*N,1:d*N) = M(1:d*N,1:d*N) + JTJ;
        
        % Bottom Left
        M(d*N + 1:d*(N+1),1:d*N) = M(d*N + 1:d*(N+1),1:d*N) + (eye(d) - J(:,:,j)*P)' * J(:,:,j);
        
        % Top Right
        M(1:d*N, d*N + 1:d*(N+1)) = M(1:d*N, d*N + 1:d*(N+1)) + J(:,:,j)' * (eye(d) - J(:,:,j)*P);
        
        % Bottom right block
        M(d*N+1:d*(N+1), d*N+1:d*(N+1)) = M(d*N+1:d*(N+1), d*N+1:d*(N+1)) + (eye(d) - J(:,:,j)*P)'*(eye(d) - J(:,:,j)*P);
    end
    
    % Bottom right block
    % M(d*N+1:d*(N+1), d*N+1:d*(N+1)) = m * eye(d);
end