function biem
clf;
[V,F] = readOBJ('plane.obj');
% [V,F] = readOBJ([data_dir() '/meshes_obj/star.obj']);
V = V(:,1:2)*1;

% x0  = [0 0; 0 2;  2 2; 2 0; ]';
% E = [1 2; 2 3; 3 4; 4 1];

E = boundary_faces(F);
% Get undeformed boundary vertices
V_bnd = unique(E(:));
x0 = V(V_bnd,:)';


% changem -- remap vertex indices
new_Idx = find(V_bnd);
[toreplace, bywhat] = ismember(E, V_bnd);
E(toreplace) = new_Idx(bywhat(toreplace));
    

% E_lines = plot([x0(1,E(:,1)); x0(1,E(:,2))], [x0(2,E(:,1)); x0(2,E(:,2))],'LineWidth',3);
axis equal
% ylim([-1 3])
hold on;

n = size(x0,2);
m = size(V,1);
N = zeros(size(x0));
y = (x0(:,E(:,2)) + x0(:,E(:,1))) / 2.0;

% Compute normals
for i = 1:n
    vdif = x0(:,E(i,2)) - x0(:,E(i,1));
    vdif = [vdif' 0];
    normal = cross(vdif, [0 0 1]);
    normal = normal / norm(normal);
    N(:,i) = normal(1:2);
    %plot([x0(1,E(i,2)); y(1,E(i,2))+normal(1)], [x0(2,E(i,2)); x0(2,E(i,2))+normal(2)])
    plot([y(1,i); y(1,i)+normal(1)], [y(2,i); y(2,i)+normal(2)])
    hold on
end

num_lines=size(y,2);
num_lines = 10;
fac=2;

xrange=[0.2 1.5];
yval=.2;
xval=0.2;
% y = [linspace(1e-4,2-1e-4,num_lines) ; repelem(0, num_lines)];
y  = [linspace(xrange(1),xrange(2),num_lines*fac); repelem(yval, num_lines*fac)];
y2 = [repelem(xval, num_lines*fac); linspace(xrange(1),xrange(2),num_lines*fac)];

% Legendre weighted points
[y_g,w_y]=lgwt(num_lines,xrange(1),xrange(2));
y_g=[y_g'; repelem(yval, num_lines)];

[y2_g,w2_y]=lgwt(num_lines,xrange(1),xrange(2));
y2_g=[repelem(xval, num_lines); y2_g'];

% y = [y y2];
% y_g = [y_g y2_g];
% w_y = [w_y; w2_y];

% Boundary integral equation
% Ref: https://appliedmath.ucmerced.edu/sites/appliedmath.ucmerced.edu/files/page/documents/introduction_to_boundary_integral_equation_methods.pdf
% Solving Lu = f where u is the densities for each boundary point
n = num_lines;
f = ones(n*fac,1);%[ones(n*fac,1); zeros(n*fac,1)];
L=zeros(n*fac*1,n*1);
inv_twopi = 1 / (2 * pi);
for i = 1:n*fac*1
    diff = y(:,i) - y_g;

    % Single layer formulation
    L(i,:) = -inv_twopi * log(vecnorm(diff)) .* w_y';
end
%mu = L \ f;
mu = pinv(L)* f;
mu = mu';

ub = zeros(n,1);
for i = 1:n
    diff = y(:,i) - y_g;
    
    % Single layer
    green = -inv_twopi * log(vecnorm(diff)) .* mu .* w_y';
    ub(i) = sum(green);
end

u = zeros(m,1);
for i = 1:m
    diff = V(i,:)' - y_g;
    
    % Single layer
    green = -inv_twopi * log(vecnorm(diff)) .* mu .* w_y';
    u(i) = sum(green);
end
max(u)
min(u)
cmap=colormap(hot(256));
% C = u ./ max(u);
u = max(u,0);
C = floor(rescale(u)*255) + 1;
C = cmap(C,:);
% C = [u u u];
p=patch('Faces',F,'Vertices',V, 'FaceVertexCData',C,'FaceColor','interp','EdgeColor','none');
hold on;
plot(y(1,:),y(2,:),'Color','m', 'LineWidth',4)
hold on;
plot(V(:,1),V(:,2),'.','Color','b');
end