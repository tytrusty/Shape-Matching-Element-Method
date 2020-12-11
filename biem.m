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
num_lines = 4;
fac=20;

% y = [linspace(1e-4,2-1e-4,num_lines) ; repelem(0, num_lines)];
y = [linspace(0.2,0.8,num_lines) ; repelem(0.3, num_lines)];
% [y_g,w_y]=lgwt(num_lines*fac,0,2);
[y_g,w_y]=lgwt(num_lines*fac,0.2,0.8);
% y_g=[y_g'; repelem(0, num_lines*fac)];
y_g=[y_g'; repelem(0.3, num_lines*fac)];
N = [repelem(0, num_lines); repelem(1, num_lines)];
n = size(N,2);

% Boundary integral equation
% Ref: https://appliedmath.ucmerced.edu/sites/appliedmath.ucmerced.edu/files/page/documents/introduction_to_boundary_integral_equation_methods.pdf
% Solving Lu = f where u is the densities for each boundary point
f = ones(n*fac,1);
%L = eye(n) * 0.5;
L=zeros(n*fac,n);
inv_twopi = 1 / (2 * pi);
for i = 1:n*fac
    % green's function
    %diff = y(:,i) - (y + [0 1e-5]');
    diff = y_g(:,i) - y;
    %diff = y - y_g(:,i);
    
    denom = vecnorm(diff).^2;
    numerator = sum(N .* diff,1);
    green = inv_twopi * (numerator ./ denom);
    %green(i)=0;
    L(i,:) = L(i,:) + green * (2/num_lines);
    
    % Single layer formulation
    L(i,:) = -inv_twopi * log(vecnorm(diff)) * (2 / num_lines);
end
mu = pinv(L) * f; %mu = L \ f;
mu = mu';

u = zeros(m,1);
for i = 1:m
    V(i,:)
    diff = V(i,:)' - y;
    denom = vecnorm(diff).^2;
    numerator = sum(N .* diff,1);
    
    % Double layer
    green = inv_twopi * (numerator ./ denom) .* mu;
    
    % Single layer
    green = -inv_twopi * log(vecnorm(diff)) .* mu;% .* w_y';
    
    u(i) = sum(green * (2 / num_lines));
end
max(u)
cmap=colormap(hot(256));
% u(isnan(u))=1;
% C = u ./ max(u);
% C = u ./ max(u);
C = floor(rescale(u)*255) + 1;
C = cmap(C,:);
% C = [u u u];
p=patch('Faces',F,'Vertices',V, 'FaceVertexCData',C,'FaceColor','interp','EdgeColor','none');
[~,I] = max(u);
hold on;
% plot(V(I,1),V(I,2),'.','MarkerSize',20,'Color','r')
% plot([y(1,:); y(1,:)+N(1,:)/5], [y(2,:); y(2,:)+N(2,:)/5])
hold on
plot(y(1,:),y(2,:), 'Color','m', 'LineWidth',4)
hold on;
plot(V(:,1),V(:,2),'.','Color','b');
end