function test_vem_dF
    n=4;
    
    % Undeformed and deformed positions
    q=sym('q',[2*n,1]);
    assume(q,'real');
    
    q0=sym('q0',[2*n,1]);
    assume(q0,'real');
    
    % 
    
%     x = 
%     BM=sym('BM',[n,1]);
%     assume(BM,'real')
    
    % Get center of massses in deformed & deformed states
    ONE=repmat([1 0; 0 1], 1, n);
    x_com=(1/n)*ONE*q;
    x0_com =(1/n)*ONE*q0;
    
    % Create P matrix
    P=sym(zeros(2, n));
    for i=1:n
        P(:,i) = q(2*i-1:2*i)-x_com;
    end
    
    % Create Q Matrix
    Q=sym(zeros(2, n));
    for i=1:n
        Q(:,i) = q0(2*i-1:2*i)-x0_com;
    end
    
%     B = Q' * inv(Q*Q');
    B = sym('B',[n,2]);
    A = P * B;
    
    X=sym('X',[2,1]);
    M = X-x0_com;
%     M(3) = M(1)^2;
%     M(4) = M(2)^2;
%     M(5) = M(1)*M(2);
    x = A*M + x_com;
    jacobian(M,X)
    AA=diff(A,q(1))
    
end