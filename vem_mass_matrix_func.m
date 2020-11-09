function M = vem_mass_matrix_func(n)
    if nargin < 1
        n=4
    end
    % n = numel(B);
    qdot=sym('qdot',[2*n,1]);
    assume(qdot,'real');
    BM=sym('BM',[n,1]);
    assume(BM,'real')

    ONE=repmat([1 0; 0 1], 1, n);
    x_com=(1/n)*ONE*qdot;
    % XCOM_q = (1/(n/2))*ONE'*ONE*qdot;

    P=sym(zeros(2, n));
    for i=1:n
        P(:,i) = qdot(2*i-1:2*i)-x_com;
    end

    C=P*BM + x_com;
    T=C'*C;
    dT_dqdot = gradient(T,qdot);
    M = jacobian(dT_dqdot,qdot)
%     M = double(subs(M, BM, B));
%     matlabFunction(M, 'vars', {BM}, 'File','mass_matrix_n','Comments','Version: 1.0')
end