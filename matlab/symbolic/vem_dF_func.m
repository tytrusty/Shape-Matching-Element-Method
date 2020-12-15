function dF = vem_dF_func(n,ne)
    if nargin < 1
        n=4;
        ne = 2;
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
    E(3:4,5:6) = eye(2);

    ONE=repmat([1 0; 0 1], 1, n);
    x_com=(1/n)*ONE*qdot;

    qE = E * qdot;
    P=sym(zeros(2, ne));
    for i=1:ne
        P(:,i) = qdot(2*i-1:2*i)-x_com;
        P(:,i) = qE(2*i-1:2*i)-x_com;
    end

    C = P * B * dM_dX;
    C1=C(:)
    dT_dqdot = jacobian(C(:),qdot)
%     M = jacobian(dT_dqdot,qdot)
%     M = double(subs(M, BM, B));
%     matlabFunction(M, 'vars', {BM}, 'File','mass_matrix_36','Comments','Version: 1.0')
end