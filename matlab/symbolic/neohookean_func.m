function neohookean_func
d=3;
F=sym('F', [d,d]);
assume(F,'real');
J=det(F);
I3=trace(F'*F)/J^(2/3);

syms C D
% neohookean
psi=C*(I3-3)+D*(J-1)^2;

% stable neohookean
% I3=det(F);
% I2=trace(F'*F);
% psi=C*(I2-3)-C*(I3-1)+D*(I3-1)^2 + (C/D)^2;

g = gradient(psi,F(:));
H = hessian(psi,F(:));
% ccode(psi)
ccode(g)
% ccode(H)
% matlabFunction(g, 'vars', {F,C,D},'File','SNH_dF','Comments','Version: 1.0')
% matlabFunction(H, 'vars', {F,C,D},'File','SNH_dF2','Comments','Version: 1.0')
end