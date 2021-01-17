function vem_simulate_nurbs_with_collision(parts, varargin)
    % Simulation parameter parsing
    p = inputParser;
    addParameter(p, 'dt', 0.01);                % timestep
    addParameter(p, 'lambda', 0.5 * 1700);    	% Lame parameter 1
    addParameter(p, 'mu', 0.5 * 15000);        	% Lame parameter 2
    addParameter(p, 'gravity', -300);          	% gravity force (direction is -z direction)
    addParameter(p, 'k_stability', 1e5);       	% stiffness for stability term
    addParameter(p, 'order', 1);              	% (1 or 2) linear or quadratic deformation
    addParameter(p, 'rho', 1);                 	% per point density (currently constant)
    addParameter(p, 'save_output', 0);        	% (0 or 1) whether to output images of simulation
    addParameter(p, 'save_obj', 0);          	% (0 or 1) whether to output obj files
    addParameter(p, 'save_resultion', 20);    	% the amount of subdivision for the output obj
    addParameter(p, 'pin_function', @(x) 1);
    addParameter(p, 'sample_interior', 1);
    addParameter(p, 'distance_cutoff', 20);
    addParameter(p, 'enable_secondary_rays', true);
    addParameter(p, 'fitting_mode', 'global');
    addParameter(p, 'plot_points', false);
    addParameter(p, 'plot_com', true);
    addParameter(p, 'com_threshold', 100000);
    addParameter(p, 'collision_ratio', 0.1);    % parameter for the collision penalty force
    addParameter(p, 'self_collision', false);
    addParameter(p, 'collision_with_other', true);
    addParameter(p, 'collision_with_plane', true);
    addParameter(p, 'collision_plane_z', -10.0);
    parse(p,varargin{:});
    config = p.Results;
    
    d = 3;  % dimension (2 or 3)
    n = numel(parts);	% number of shapes
    
    % The number of elements in the monomial basis.
    k = basis_size(d, config.order);
    
    % Read in NURBs 
    fig=figure(1);
    clf;
    parts=nurbs_plot(parts);

    % Assembles global generalized coordinates
    [J, hires_J, q, E, x0] = nurbs_assemble_coords(parts);

    % Generating centers of mass. Temporary method!
    [x0_coms, com_cluster, com_map] = generate_com(parts, x0, E, ...
        config.com_threshold, n);
    
    if config.plot_com
        com_plt = plot3(x0_coms(1,:),x0_coms(2,:),x0_coms(3,:), ...
                        '.','Color','g','MarkerSize',20);
        hold on;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%
    
    % Initial deformed positions and velocities
    x = x0;
    qdot=zeros(size(q));
    
    % Setup pinned vertices constraint matrix
    pin_I = config.pin_function(x0);
    P = fixed_point_constraint_matrix(x0',sort(pin_I)');
    
    % Plotting pinned vertices.
    X_plot=plot3(x(1,pin_I),x(2,pin_I),x(3,pin_I),'.','Color','red','MarkerSize',20);
    hold on;
    
    % Sampling points used to compute energies.
    if config.sample_interior
        [V, vol] = raycast_quadrature(parts, [9 9], 5);
    else
    	V=x0;
        vol=ones(size(V,2),1);
    end
    m = size(V,2);  % number of quadrature points

    if config.plot_points
        V_plot=plot3(V(1,:),V(2,:),V(3,:),'.','Color','m','MarkerSize',20);
    end

    % Lame parameters concatenated.
    params = [config.mu * 0.5, config.lambda * 0.5];
    params = repmat(params,size(V,2),1);
        
    % Gravity force vector.
  	f_gravity = repmat([0 0 config.gravity], size(x0,2),1)';
    f_gravity = config.dt*P*f_gravity(:);
    
    % Undeformed Center of mass
    x0_com = mean(x0,2);
        
    % Shape Matrices
    L = compute_shape_matrices(x0, x0_coms, com_map, E, ...
        com_cluster, config.order, config.fitting_mode);
    
    % Compute Shape weights
    [w, w_I] = nurbs_blending_weights(parts, V', config.distance_cutoff, ...
        'Enable_Secondary_Rays', config.enable_secondary_rays);
    [w0, w0_I] = nurbs_blending_weights(parts, x0', config.distance_cutoff, ...
        'Enable_Secondary_Rays', config.enable_secondary_rays);
                             
    % Build Monomial bases for all quadrature points
    [Y,Y_S] = vem_dx_dc(V, x0_coms, w, w_I, com_map, config.order, k);
    [Y0,Y0_S] = vem_dx_dc(x0, x0_coms, w0, w0_I, com_map, config.order, k);
    
    % Fixed x values.
    x_fixed = zeros(size(x0));
    for i = 1:numel(pin_I)
        x_fixed(:,pin_I(i))=x0(:,pin_I(i));
    end
    
    % Applying fixed point constraints to NURBS jacobian.
    J = P * J;
    
    % Computing each gradient of deformation gradient with respect to
    % projection operator (c are polynomial coefficients)
    [dF_dc, dF_dc_S] = vem_dF_dc(V, x0_coms, w, w_I, com_map, config.order, k);

    % Compute mass matrices
    ME = vem_error_matrix(Y0, Y0_S, L, d);
    M = vem_mass_matrix(Y, Y_S, L, config.rho .* vol);
    M = (M + config.k_stability*ME); % sparse?
    % Save & load these matrices for large models to save time.
    % save('saveM.mat','M');
    % save('saveME.mat','ME');
    % M = matfile('saveM.mat').M;
    % ME = matfile('saveME.mat').ME;
    
    % Triangulate nurbs patch, for collision detection
    [verts,faces] = triangulate_iges(parts);
    collision_ratio = config.collision_ratio;
    face_normal = normals(verts, faces);
    
    if config.collision_with_other
      collide_iges_file = 'rounded_cube.iges';
      collide_parts = nurbs_from_iges(collide_iges_file);
      [collide_verts, collide_faces] = triangulate_iges(collide_parts);
      collide_verts = collide_verts * (max(verts(:,1))-min(verts(:,1))) / (max(collide_verts(:,1))-min(collide_verts(:,1)))  + [0 0 max(verts(:,3))];
      t_collide = tsurf(collide_faces, collide_verts);
      move_ratio = (max(collide_verts(:,1))-min(collide_verts(:,1)));
      set(fig, 'KeyPressFcn', @keypress)
    end
    
    % draw the plane if needed
    if config.collision_with_plane
      if config.collision_plane_z > min(verts(:,3))
        config.collision_plane_z = min(verts(:,3)) - 20.0;
      end
      hold on
      x_range = 2 * (max(verts(:,1))-min(verts(:,1)));
      y_range = 2 * (max(verts(:,2))-min(verts(:,2)));
      p = patch('XData', x_range*[-1 -1 1 1],'YData', y_range*[-1 1 1 -1], ...
                'ZData', config.collision_plane_z * [1 1 1 1], ...
                'FaceColor',[0.5 0.5 0.5], 'FaceAlpha', 0.5, 'FaceLighting', 'gouraud');
      axis equal
      hold off
    end

    ii=1;
    for t=0:config.dt:30
        tic
        
        % Preparing input for stiffness matrix mex function.
        % TODO: don't form this vector this way :)
        b = [];
        for i=1:numel(E)
            b = [b (x(:,E{i}))];
        end
        b = b(:);

        % Solve for polynomial coefficients (projection operators).
        c = L * b;
        p = c(end-d+1:end);
        x_com = x0_com + p;   
        
        % Stiffness matrix (mex function)
        K = -vem3dmesh_neohookean_dq2(c, vol, params, dF_dc, w_I, k, n, ...
                                      size(x0_coms,2));
        K = L' * K * L;
        
        % Force vector
        dV_dq = zeros(d*(k*n + size(x0_coms,2)),1);

        % Computing force dV/dq for each point.
        % TODO: move this to C++ :)
        for i = 1:m
            % Deformation Gradient
            F = dF_dc{i} * dF_dc_S{i} * c;
            F = reshape(F,d,d);
            
            V(:,i) = Y{i} * Y_S{i} * c;
            
            % Force vector
            dV_dF = neohookean_tet_dF(F, params(i,1), params(i,2));
            dV_dq = dV_dq +  dF_dc_S{i}' * dF_dc{i}' * dV_dF * vol(i);
        end
        dV_dq = L' * dV_dq;
        
        % Error correction force
        f_error = - 2 * ME * x(:);
        
        f_error = config.k_stability*(config.dt * P * f_error(:));
       
        % Force from potential energy.
        f_internal = -config.dt*P*dV_dq;
        
        % Force for collision.
        f_collision = zeros(size(verts,1)*3,1);
        % detect the collision with another mesh
        if config.collision_with_other
          [IF] = intersect_other(verts,faces,collide_verts,collide_faces);
          for i = 1:size(IF,1)
            fid = IF(i, 1); % colliding face id of the origin mesh
            c_fid = IF(i, 2); % colliding face id of the colliding mesh
            f_collision = check_collision_between_faces(verts,faces,collide_verts,collide_faces,fid,c_fid,face_normal,collision_ratio,f_collision);
          end
        end
        % self-collision detection
        if config.self_collision
          [~,~,IF] = selfintersect(verts,faces,'DetectOnly',true);
          for i = 1:size(IF,1)
            fid = IF(i, 1); % colliding face id of the origin mesh
            c_fid = IF(i, 2); % colliding face id of the colliding mesh
            if size(intersect(faces(fid,:), faces(c_fid,:)), 2) == 0 
              continue; % skip the face pairs if they share a vertex
            end
            f_collision = check_collision_between_faces(verts,faces,verts,faces,fid,c_fid,face_normal,collision_ratio,f_collision);
            f_collision = check_collision_between_faces(verts,faces,verts,faces,c_fid,fid,face_normal,collision_ratio,f_collision);
          end
        end
        % collision with the plane
        if config.collision_with_plane
          [collide_plane_vid, ~] = find(verts(:,3) < config.collision_plane_z);
          f_collision_plane = zeros(size(verts,1), 3);
          f_collision_plane(collide_plane_vid, 3) = 10 * collision_ratio * ...
                                                    (config.collision_plane_z*ones(size(collide_plane_vid,1), 1) - verts(collide_plane_vid,3));
          f_collision_plane = reshape(f_collision_plane', [], 1);                           
          f_collision = f_collision + f_collision_plane;
        end
%         B = barycenter(verts, faces);
%         for i = 1:size(IF,1)
%           fid = IF(i, 1); % colliding face id of the origin mesh
%           c_fid = IF(i, 2); % colliding face id of the colliding mesh
%            if size(intersect(faces(fid,:), faces(c_fid,:)), 2) == 0 
%              continue;
%            end
%           hold on
%           plot3(B([fid;c_fid],1),B([fid;c_fid],2),B([fid;c_fid],3));
%         end
%         hold off
        f_collision = hires_J' * f_collision;
        
        % Computing linearly-implicit velocity update
        lhs = J' * (P*(M - config.dt*config.dt*K)*P') * J;
        rhs = J' * (P*M*P'*J*qdot + f_internal + f_gravity + f_error) + f_collision;
        qdot = lhs \ rhs;

        % Update position
        q = q + config.dt*qdot;
        x = reshape(P'*J*q,3,[]) + x_fixed;
        
        % Update high-resolution triangulation
        verts = reshape(hires_J * q, 3, [])';
        face_normal = normals(verts, faces);
        
        % Update NURBs plots
        x_idx=0;
        for i=1:numel(parts)
            x_sz = size(parts{i}.x0,2);
            xi = x(:,x_idx+1:x_idx+x_sz);
            parts{i}.plt.Vertices =xi';
            x_idx = x_idx+x_sz;
        end
        
        if config.plot_com
            x_coms = c(d*k*n + 1:end); % extract centers of mass
            com_plt.XData = x_coms(1:d:end);
            com_plt.YData = x_coms(2:d:end);
            com_plt.ZData = x_coms(3:d:end);
        end
        
        if config.plot_points
            V_plot.XData = V(1,:);
            V_plot.YData = V(2,:);
            V_plot.ZData = V(3,:);
        end
        drawnow
        
        if config.save_obj
            obj_fn = "output/obj/part_" + int2str(ii) + ".obj";
            nurbs_write_obj(q,parts,obj_fn,ii);
        end
        
        if config.save_output
            fn=sprintf('output/img/cutoff_ten_%03d.png',ii);
            saveas(fig,fn);
        end
        ii=ii+1
        toc
    end
    
    % Callback to process keypress events
    function keypress(~, evnt)
      switch lower(evnt.Key)  
          case 'leftarrow'
             collide_verts = collide_verts - 0.1 * move_ratio * repmat([1 0 0], size(collide_verts,1), 1);
          case 'rightarrow'
             collide_verts = collide_verts + 0.1 * move_ratio * repmat([1 0 0], size(collide_verts,1), 1);
          case 'downarrow'
             collide_verts = collide_verts - 0.1 * move_ratio * repmat([0 1 0], size(collide_verts,1), 1);
          case 'uparrow'
             collide_verts = collide_verts + 0.1 * move_ratio * repmat([0 1 0], size(collide_verts,1), 1);
          case 115  % s: move down
             collide_verts = collide_verts - 0.1 * move_ratio * repmat([0 0 1], size(collide_verts,1), 1);
          case 119  % w: move up
             collide_verts = collide_verts + 0.1 * move_ratio * repmat([0 0 1], size(collide_verts,1), 1);
          case 97   % a: increase the size of the sphere
             collide_verts = collide_verts * 1.1;
          case 100  % d: decrease the size of the sphere
             collide_verts = collide_verts * 0.9;
          otherwise
              return
      end
      % Always do a redraw
      t_collide.Vertices = collide_verts;
      drawnow
    end
end

function [verts,faces] = triangulate_iges(parts)
  faces=[];
  verts=[];
  for ii=1:numel(parts)
      if isfield(parts{ii}, 'hires_T')
          F = parts{ii}.hires_T;
          V = parts{ii}.hires_x0';
      else
          F = parts{ii}.T;
          V = parts{ii}.x0';
      end
      F = F + size(verts,1);
      faces=[faces; F];
      verts=[verts; V];
  end
  % check and remove degenerate faces
  dblA = doublearea(verts,faces);
  deg_fid = find(dblA <= 0);
  faces(deg_fid, :) = [];
end

function f_collision = check_collision_between_faces(verts,faces,collide_verts,collide_faces,fid,c_fid,face_normal,collision_ratio,f_collision)
  % check vertex-face collisions in the colliding faces
  for j = 1:3  
    cv = collide_verts(collide_faces(c_fid,j),:);
    [dist, cp] = pointTriangleDistance(...
                            [verts(faces(fid,1),:); verts(faces(fid,2),:); verts(faces(fid,3),:)], ...
                            cv);
    if dot(cv-cp, face_normal(fid, :)) > 0
      f_collision(3*faces(fid,1)-2:3*faces(fid,1), 1) = f_collision(3*faces(fid,1)-2:3*faces(fid,1), 1) - 1/3 * collision_ratio * dist * face_normal(fid,:)';
      f_collision(3*faces(fid,2)-2:3*faces(fid,2), 1) = f_collision(3*faces(fid,2)-2:3*faces(fid,2), 1) - 1/3 * collision_ratio * dist * face_normal(fid,:)';
      f_collision(3*faces(fid,3)-2:3*faces(fid,3), 1) = f_collision(3*faces(fid,3)-2:3*faces(fid,3), 1) - 1/3 * collision_ratio * dist * face_normal(fid,:)';
      disp('face collide!!');
    end
  end
end


