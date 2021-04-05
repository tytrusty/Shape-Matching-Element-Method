function neohookean_func
d=3;
F=sym('F', [d,d]);
assume(F,'real');
J=det(F);
I3=trace(F'*F)/J^(2/3);

syms C D
% neohookean
psi=C*(I3-3)+D*(J-1)^2;

% % stable neohookean
% I3=det(F);
% I2=trace(F'*F);
% psi=C*(I2-3)-2*C*(I3-1)+D*(I3-1)^2;

% SNH
syms Mu La
J = det(F);
IC = trace(F'*F);
a = 1.0 + (Mu/La) - (Mu / (4*La));
psi = Mu*(IC-3) + La*(J-a)^2 - Mu*log(IC+1);

g = gradient(psi,F(:));
H = hessian(psi,F(:));
% ccode(psi)
ccode(g)
ccode(H)
matlabFunction(g, 'vars', {F,Mu,La},'File','SNH_dF','Comments','Version: 1.0')

% matlabFunction(g, 'vars', {F,C,D},'File','SNH_dF','Comments','Version: 1.0')
% matlabFunction(H, 'vars', {F,C,D},'File','SNH_dF2','Comments','Version: 1.0')
end