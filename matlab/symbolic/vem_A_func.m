function A = vem_A_func(n,ne)
    if nargin < 1
        n=4;
        ne = 2;
    end
    % n = numel(B);
    qdot=sym('qdot',[2*n,1]);
    assume(qdot,'real');
    BM=sym('BM',[ne,1]);
    assume(BM,'real')
    
    E = zeros(2*ne,2*n);
    E(1:2,3:4) = eye(2);
    E(3:4,5:6) = eye(2);

    ONE=repmat([1 0; 0 1], 1, n);
    x_com=(1/n)*ONE*qdot;
    x_com = sym('x_com',[2,1]);
    
    qE = E * qdot;
     P=sym(zeros(2, ne));
    for i=1:ne
        %P(:,i) = qdot(2*i-1:2*i)-x_com;
        %P(:,i) = qE(2*i-1:2*i) - x_com;
        P(:,i) = qE(2*i-1:2*i)- [qdot(end-1) qdot(end)]' - x_com;
    end
    
%     C=([qdot(3) qdot(4)]' -  (P*BM + x_com)); % Error matrix term
    C=(P*BM + x_com);                     % Mass matrix term
    dT_dqdot = jacobian(C,qdot)
%     M = jacobian(dT_dqdot,qdot)
%     M = double(subs(M, BM, B));
%     matlabFunction(M, 'vars', {BM}, 'File','mass_matrix_36','Comments','Version: 1.0')
end