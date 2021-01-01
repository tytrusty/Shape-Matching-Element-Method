function vem_test
    if nargin < 1
        n=4;
        ne = 1;
        k=5;
    end
    % n = numel(B);
    qdot=sym('qdot',[2*n,1]);
    assume(qdot,'real');
    
    
    B=sym('B',[ne,k]);
    assume(B,'real')
    
    dM_dX=sym('dM_dX',[k,2]);
    assume(dM_dX,'real');    
    
    E = zeros(2*ne,2*n);
    E(1:2,3:4) = eye(2);
%     E(3:4,5:6) = eye(2);

    x_com = sym('x_com',[2,1]);
    qE = E * qdot;
    P=sym(zeros(2, ne));
    for i=1:ne
%         P(:,i) = qdot(2*i-1:2*i)-x_com;
        P(:,i) =   qE(2*i-1:2*i)- [qdot(end-1) qdot(end)]';
    end

    C = P * B * dM_dX;
    C1=C(:);
    dT_dqdot = jacobian(C(:),qdot)
end