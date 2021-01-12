%cantileaver simulation
[V,T,F] = readMESH('vem_beam_tet2.mesh');


%setup simulation variables
   
%tet volumes
vol = volume(V,T);

%triangle gradients 
dphidX = linear_tetmesh_dphi_dX(V,T);

%Mass Matrix
rho =  100.0*ones(size(T,1),1);
M = linear_tetmesh_mass_matrix(V,T, rho, vol);

%material properties 
YM = 5e5; %in Pascals
pr =  0.45;
[lambda, mu] = emu_to_lame(YM*ones(size(T,1),1), pr*ones(size(T,1),1));

%boundary conditions are fun
min_I = find(V(:,1) == min(V(:,1)));
P = fixed_point_constraint_matrix(V, sort(min_I));

%external forces
gravity = [0 -10 0]'; % We are on krypton 
gravity = P*M*repmat(gravity, size(V,1),1);

%project mass matrix down
M = P*M*P';

%warm up the display
[~,p] = nice_plot(V,F);

%simulation setup and loop
energy_func = @(a,b, c, d, e) linear_tetmesh_neohookean_q(a,b,c,d,e,[0.5*mu, 0.5*lambda]);
gradient_func = @(a,b, c, d, e) linear_tetmesh_neohookean_dq(a,b,c,d,e,[0.5*mu, 0.5*lambda]);
hessian_func = @(a,b, c, d, e) linear_tetmesh_neohookean_dq2(a,b,c,d,e,[0.5*mu, 0.5*lambda], 'fixed');

dt = 0.1;

qt = reshape(V', 3*size(V,1),1);

%setup all thw variables I need using bc projection matrix
b = qt - P'*P*qt;
qt = P*qt;
vt = 0*qt;


disp('Press SPACEBAR to CONTINUE');
pause 

h = text(0.5*(max(V(:,1)) + min(V(:,1))), 0.9*max(V(:,2)), num2str(0));
h.FontSize = 28;
    
for ti=1:500
    
     h.String = num2str(ti);

     g =  -M*vt + ...
           dt*P*gradient_func(V,T, P'*qt+b, dphidX, vol) + ...
           - dt*gravity;
       
     H = M + dt*dt*P*hessian_func(V,T, P'*qt+b, dphidX, vol)*P';

     vt = -H\g; 
     qt = qt + dt*vt;

     p.Vertices = reshape(P'*qt + b, 3, size(V,1))';
     drawnow;
end
