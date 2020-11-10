function neohookean_func
F=sym('F', [2,2]);
assume(F,'real');
J=det(F);
I3=trace(F'*F)/J^(2/3);

syms C D
psi=C*(I3-3)+D*(J-1)^2;

g = gradient(psi,F(:))
H = hessian(psi,F(:));
% matlabFunction(g, 'vars', {F,C,D},'File','neohookean_dF','Comments','Version: 1.0')
matlabFunction(H, 'vars', {F,C,D},'File','neohookean_dF2','Comments','Version: 1.0')
end