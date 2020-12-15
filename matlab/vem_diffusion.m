%just doing some shape matching plus figuring out how to flatten
%symmetric tensors

global v_d;
global V;
global T;
global areas;
global dphidX;

v_d = [];

%load a mesh
[V,T] = readOBJ([data_dir() '/meshes_obj/star.obj']);
 V(1,:) = [];
 T = T - 1;
 
 areas = triangle_area(V,T);
 dphidX = linear_tri2dmesh_dphi_dX(V(:,1:2),T);
 
 G  = grad(V,T);
 dblA = doublearea(V,T);
 L = G'*repdiag(diag(sparse(dblA)/2),size(V,2))*G;
    
global F;
global F_local

F = boundary_faces(T);

%F = F(1,:);
%build the virtual shape functions

%1. compute center of mass 
global X_c;
X_c = sum(V,1)./size(V,1);

v = V;

%1b weighting functions for weight computation 

%2. shape matching to each edge
%X = nx3 try to stick with igl format as long as possible because I don't
%know
Xe1 = V(F(:,1),:);
Xe2 = V(F(:,2),:);

%this should be extracted to a utility function (block matrix to sparse
%matrix)
Dtmp = full([kroneye(reshape((Xe1 - X_c)', 1, numel(Xe1)), 3); kroneye(reshape((Xe2 - X_c)', 1, numel(Xe2)), 3)]);

jcol = repmat(1:(size(F,1)*9), 6, 1);
jcol = jcol(:);

icol = repmat(reshape(1:(size(F,1))*6, 6, size(F,1)), 9,1);
icol = icol(:);

D = sparse(icol, jcol, Dtmp(:));

%minimize difference between monomial and points of each boundary edge
ve1 = v(F(:,1), :) - X_c;
ve2 = v(F(:,2), :) - X_c;

%projection for 2d
global P;
P = speye(9*size(Xe1,1),9*size(Xe1,1));
P(:, 9:9:end) =  [];
P(:, 8:8:end) =  [];
P(:, 7:7:end) =  [];
P(:, 6:6:end) =  [];
P(:, 3:5:end) =  [];

%s are monomial coefficients
s = (P'*D'*D*P)\(P'*D')*(reshape([ve1';ve2'], 2*numel(ve1),1));
s = P*s;

%ugly but just storing some global variables here for doing this solve
%interactively 
global A_mat;
global B_mat;
A_mat = (P'*D'*D*P)
B_mat = (P'*D');

%plot some shapes 
v_d = (V - X_c)*reshape(s(28:36),3,3)' + X_c;

%try Ty's distance based function 

%for each point, find closes point on an edge 

%edge vectors 
dE = (V(F(:,2),:) - V(F(:,1),:));
norm2_dE = sum(dE.*dE, 2);

dVx = V(:,1)' - V((F(:,1)),1);
dVy = V(:,2)' - V((F(:,1)),2);
dVz = V(:,3)' - V((F(:,1)),3);

dEdotdV = (dE(:,1).*dVx +  dE(:,2).*dVy +dE(:,3).*dVz)./norm2_dE;

%bound constraints
dEdotdV = min(max(dEdotdV, 0),1);

point_index = 10;
VE = dEdotdV(:,point_index).*dE + V(F(:,1),:);

global a;
a = zeros(size(V,1), size(F,1));
w = zeros(size(V,1), size(F,1));

%solve for one function to see what it looks
%plot th;e shape functions

%pick 2 edges and build shape function 

dist_edge = sqrt(((dEdotdV.*dE(:,1) + V(F(:,1),1)) - V(:,1)').^2 + ((dEdotdV.*dE(:,2) + V(F(:,1),2)) - V(:,2)').^2 + ((dEdotdV.*dE(:,3) + V(F(:,1),3)) - V(:,3)').^2)';

%diffusion weighting function for each edge
for jj=1:size(F,1)
    
    bc = zeros(size(V,1),1);
    bc(F(jj,1)) = 1;
    bc(F(jj,2)) = 1;
    bc = L*bc; 
    
    L_mod = L;
    L_mod(:, unique(F(:))) = [];
    L_mod(unique(F(:)),:) = [];
    bc(unique(F(:))) = [];
    w(F(jj,1), jj) = 1;
    w(F(jj,2), jj) = 1;
    
    unique(F(:)) %all nodes on edges, don't put anything here
    w(setdiff(1:size(V,1), unique(F(:))), jj) = -L_mod\bc;
end

max(w,[],'all')
%w = exp(-dist_edge.*dist_edge); %basis 


for ii = 1:size(a,1)
    
    W = diag(1./w(ii,:));
    W(isinf(W)) = 1e8;
    
    a(ii,:) = quadprog(W, zeros(size(F,1),1), [], [],...
                        ones(1,size(F,1)), 1,...
                       zeros(size(F,1),1), ones(size(F,1),1));
    
end

%use diffusion weighting to ocmpute each point
% for ii=1:size(a,1)
%     VE = dEdotdV(:,ii).*dE + V(F(:,1),:);
%     a(ii,:) = quadprog(diag(w(ii,:).^2), zeros(size(F,1),1), [diag(w(ii,:)); -diag(w(ii,:))], [ones(size(a,2),1) zeros(size(a,2),1)],...
%                        [w(ii,:);VE(:,1)'; VE(:,2)'; VE(:,3)'] , [1; V(ii,:)'],...
%                       [], []);
%     
%     %a(ii,:) = linprog(ones(size(F,1),1), [], [],...
%     %                   [ones(1,size(F,1));VE(:,1)'; VE(:,2)'; VE(:,3)'] , [1; V(ii,:)'],...
%     %                   zeros(size(F,1),1), ones(size(F,1),1));
%                    
% end
 

%final weighting functions
%a = w;

%deform the mesh to see what happens

%make a little UI to play with
f = figure;
scale = 3;
hold on
global h;
h = tsurf(T, v_d);
x_lower = min(V,[],1);
x_upper = max(V,[],1);
x_c = [0.5*(x_lower(1) + x_upper(1)), 0.5*(x_lower(2) + x_upper(2))];
f.CurrentAxes.XLim = [scale*(x_lower(1)-x_c(1))+x_c(1), scale*(x_upper(1)-x_c(1))+x_c(1)];
f.CurrentAxes.YLim = [scale*(x_lower(2)-x_c(2))+x_c(2), scale*(x_upper(2)-x_c(2))+x_c(2)];

%set callback functions
f.WindowButtonDownFcn = @button_down_callback;
f.WindowButtonUpFcn = @button_up_callback;
%f.WindowButtonMotionFcn = @button_motion_callback;
%addlistener(gca,'WindowButtonMotionEvent',...
 %                            @(hSource,eventData) button_motion_callback(app, event, V, F));
                         
hold off

%functions and variables to enable interaction
global selected_id;
global x_selected;
selected_id = [];
x_selected = [0 0 0];


%v is vertices of edge positions
function cost = deformation_cost(v)
    
    global A_mat;
    global B_mat;
    global F;
    global X_c;
    global P;
    global V;
    global T;
    global a;
    global areas;
    global dphidX;
    
    YM = 5e6; %in Pascals
    pr =  0.45;
    [lambda, mu] = emu_to_lame(YM*ones(size(T,1),1), pr*ones(size(T,1),1));
    
    v_d = 0*V;
    v_d(F(:,1), :) = reshape(v,3, numel(v)/3)';
  
    ve1 = v_d(F(:,1), :) - X_c;
    ve2 = v_d(F(:,2), :) - X_c;

    %recompute shape matching stuff 
    %s are monomial coefficients
    s = A_mat\B_mat*(reshape([ve1';ve2'], 2*numel(ve1),1));
    s = P*s;

    
    v_tmp = 0*v_d;
    %plot some shapes 
    for ii = 1:9:numel(s)
        v_tmp = v_tmp + a(:, (ii-1)/9 +1).*((V - X_c)*reshape(s(ii:(ii+8)),3,3)');
    end
    
    v_d = v_tmp + X_c; 
    
    %compute energy of triangle mesh
    cost = linear_tri2dmesh_neohookean_q(V(:,1:2),T, igl2bart(v_d(:,1:2)), dphidX, areas, [0.5*mu, 0.5*lambda]);
    
end

%select a point on mouse click
function button_down_callback(app, event)
    global selected_id;
    global x_selected;
    global v_d;
    global F;
    
     if isMultipleCall()
        disp('blah')
        return
     end
    
    if ~isempty(selected_id)
        return;
    end
    
    x_selected = event.Source.CurrentAxes.CurrentPoint(1,:); %for 2D view
   
    %find nearest
    [~, selected_id] = min((v_d(:,1) -x_selected(1)').^2 + (v_d(:,2) - x_selected(2)').^2);
    
    if ~isempty(selected_id)
       if sum(selected_id == F(:)) == 0
           selected_id = [];
       end
    end
    
    f = gcf();
    set(f,'WindowButtonDownFcn',[])
    set(f,'WindowButtonMotionFcn',@button_motion_callback)
    set(f,'WindowButtonUpFcn',@button_up_callback)
    
    selected_id
end
%right way to implement all this is to take in an object that has
%approriate properties, that's not toooooo hard
%for now, do it all hacky and gross using global variables
%ask forgiveness later


function button_up_callback(app, event, handles)
    global selected_id;
    global A_mat;
    global B_mat;
    global v_d;
    global F;
    global X_c;
    global P;
    global V;
    global h;
    global a;
    global V;
    global T;
    global areas;
    global dphidX;
    
    f = gcf;
    set(f,'WindowButtonUpFcn',[])
    set(f,'WindowButtonMotionFcn',[])
    disp('Up');
    
    tmp_v = v_d;
    J = zeros(1, size(v_d,1));
    J(selected_id) = 1;
    J = J(F(:,1));
    J = kron(J, eye(3,3));
    b  = tmp_v(selected_id, :)';
    
    selected_id = [];
    
    options = optimoptions('fmincon');
    %options.Algorithm = 'trust-region-reflective';
    options.Display = 'iter';
    %options.SpecifyObjectiveGradient = true;
    options.MaxIterations = 100;
    options.MaxFunctionEvaluations = 1e6;
    options.FiniteDifferenceType= 'central';
    options.HessianApproximation = 'bfgs';
    %options.UseParallel = true;
    
    y = fmincon(@deformation_cost, igl2bart(tmp_v(F(:,1),:)), [], [], J, b, [], [], [], options);
    
    tmp_v(F(:,1),:) = bart2igl(y,3);
    
    ve1 = tmp_v(F(:,1),:) - X_c;
    ve2 = tmp_v(F(:,2),:) - X_c;
 
     %recompute shape matching stuff 
     %s are monomial coefficients
     s = A_mat\B_mat*(reshape([ve1';ve2'], 2*numel(ve1),1));
     s = P*s;

     v_tmp = 0*v_d;
     %plot some shapes 
     for ii = 1:9:numel(s)
         v_tmp = v_tmp + a(:, (ii-1)/9 +1).*((V - X_c)*reshape(s(ii:(ii+8)),3,3)');
     end
    
     v_d = v_tmp + X_c; 
    
      YM = 5e6; %in Pascals
    pr =  0.45;
    [lambda, mu] = emu_to_lame(YM*ones(size(T,1),1), pr*ones(size(T,1),1));
    cost = linear_tri2dmesh_neohookean_q(V(:,1:2),T, igl2bart(v_d(:,1:2)), dphidX, areas, [0.5*mu, 0.5*lambda])
    
    h.Vertices = v_d;
    drawnow limitrate nocallbacks
        
    
    
    
    set(f,'WindowButtonDownFcn',@button_down_callback);
end

function button_motion_callback(app, event)
    global selected_id;
    global x_selected;
    global v_d;
    global h;
     
    if isMultipleCall()
        disp('blah')
        return
     end
    
    if ~isempty(selected_id)
        x_selected_new = event.Source.CurrentAxes.CurrentPoint(1,:); %for 2D view
        dx  = x_selected_new - x_selected;
        v_d(selected_id,:) = v_d(selected_id,:) +dx;
        
        %move selected vertex
        x_selected = x_selected_new; 
        h.Vertices = v_d;
        
        drawnow limitrate nocallbacks
    end
   
    f = gcf;
    set(f,'WindowButtonUpFcn',@button_up_callback)
end


function flag=isMultipleCall()
s = dbstack();
% s(1) corresponds to isMultipleCall
if numel(s)<=2, 
    flag=false; 
    return; 
end
% compare all functions on stack to name of caller
count = sum(strcmp(s(2).name,{s(:).name}));
% is caller re-entrant?
if count>1, 
    flag=true; 
else
    flag=false;
end
end
%function to handle draffiing