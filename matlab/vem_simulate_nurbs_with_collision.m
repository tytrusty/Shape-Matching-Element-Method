function vem_simulate_nurbs_with_collision_new(parts, varargin)
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
    addParameter(p, 'fitting_mode', 'hierarchical');
    addParameter(p, 'plot_points', false);
    addParameter(p, 'plot_com', true);
    addParameter(p, 'collision_ratio', 0.1);             % parameter for the collision penalty force
    addParameter(p, 'self_collision', false);            % enable collision detection of the self-collision
    addParameter(p, 'collision_with_other', false);      % enable collision detection with another mesh
    addParameter(p, 'collision_with_other_sim', false);  % enable the collision response of the other mesh
    addParameter(p, 'collision_other_position', []);     % the initial position of the other mesh
    addParameter(p, 'collision_with_plane', false);      % enable collision detection with the plane
    addParameter(p, 'collision_plane_z', -10.0);         % position of the plane
    addParameter(p, 'collision_with_sphere', false);     % enable collision detection with multiple spheres
    addParameter(p, 'collision_sphere_c', []);           % center of the sphere
    addParameter(p, 'collision_sphere_r', 0.05);         % radius of the sphere
    addParameter(p, 'collision_sphere_rho', 2e3);        % density of the sphere
    addParameter(p, 'collision_sphere_initial_veloity', [0 0 0]);        % density of the sphere
    addParameter(p, 'initial_velocity', [0 0 0]);        % initial velocity of the nurbs model
    addParameter(p, 'x_samples', 5);
    addParameter(p, 'y_samples', 9);
    addParameter(p, 'z_samples', 9);
    addParameter(p, 'f_external', [0 0 0]);
    addParameter(p, 'f_external_time', 1000);     % lasting time of the external force
    addParameter(p, 'save_obj_path', 'output/obj/');
    
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
    [J, hires_J, q, E, x0, S] = nurbs_assemble_coords(parts);
    
    % Initial deformed positions and velocities
    x = x0;
%     qdot=zeros(size(q));
    qdot = reshape(repmat(config.initial_velocity, size(q,1)/3, 1)', [], 1);
    
    % Setup pinned vertices constraint matrix
    pin_I = config.pin_function(x0);
    P = fixed_point_constraint_matrix(x0',sort(pin_I)');
    
    % Plotting pinned vertices.
    X_plot=plot3(x(1,pin_I),x(2,pin_I),x(3,pin_I),'.','Color','red','MarkerSize',20);
    hold on;
    
    % Sampling points used to compute energies.
    if config.sample_interior
        yz_samples = [config.y_samples config.z_samples];
        [V, vol] = raycast_quadrature(parts, yz_samples, config.x_samples);
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

    % Compute Shape weights
    [w, w_I] = nurbs_blending_weights(parts, V', config.distance_cutoff, ...
        'Enable_Secondary_Rays', config.enable_secondary_rays);
    [w0, w0_I] = nurbs_blending_weights(parts, x0', config.distance_cutoff, ...
        'Enable_Secondary_Rays', config.enable_secondary_rays);

    % Generate centers of mass.
    [x0_coms, com_cluster, com_map] = generate_com(x0, E, w, n, parts);
    if config.plot_com
        com_plt = plot3(x0_coms(1,:),x0_coms(2,:),x0_coms(3,:), ...
                        '.','Color','g','MarkerSize',20);
        hold on;
    end
    
    % Shape Matrices
    L = compute_shape_matrices(x0, x0_coms, com_map, E, ...
        com_cluster, config.order, config.fitting_mode, S);

    % Build Monomial bases for all quadrature points
    [Y,Y_S,com_I] = vem_dx_dc(V, x0_coms, w, w_I, com_map, config.order, k);
    [Y0,Y0_S,C0_I] = vem_dx_dc(x0, x0_coms, w0, w0_I, com_map, config.order, k);
    
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
    
     % Gravity force vector.
    dg_dc = vem_ext_force([0 0 config.gravity]', config.rho*vol, Y, Y_S);
    f_gravity = config.dt*P*(L' * dg_dc);
    
    % Optional external force vector
    dext_dc = vem_ext_force(config.f_external', config.rho*vol, Y, Y_S);
    f_external = config.dt*P*(L' * dext_dc); 

    % Compute mass matrices
    ME = vem_error_matrix(Y0, L, w0_I, C0_I, d, k, n);
    M = vem_mass_matrix(Y, L, config.rho .* vol, w_I, com_I, d, k, n);

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
    %writeOBJ('./test.obj', verts, faces, [], zeros(size(faces,1),3), face_normal, zeros(size(faces,1),3));
    
    if config.collision_with_other
      collide_iges_file = 'puft_full_simpler.iges';
      collide_parts = nurbs_from_iges(collide_iges_file);
      [collide_verts, collide_faces] = triangulate_iges(collide_parts);
%       [collide_verts, collide_faces] = readOBJ('./models/fandisk.obj');
%       collide_verts = 5 * collide_verts * (max(verts(:,1))-min(verts(:,1))) / (max(collide_verts(:,1))-min(collide_verts(:,1)));
      if size(config.collision_other_position, 1) == 0 
        collide_verts = collide_verts + [0 0 max(verts(:,3)) + 10.0];
      else
        tmp = [-1  0  0; 0 -1  0; 0  0  1];
        collide_verts = (tmp * collide_verts')';
        collide_verts = collide_verts + config.collision_other_position;
      end
      t_collide = tsurf(collide_faces, collide_verts);
      move_ratio = (max(collide_verts(:,1))-min(collide_verts(:,1)));
      set(fig, 'KeyPressFcn', @keypress)
      collide_other_v = [0; 0; -10];
      collide_other_m = 1.0;
    end
    
    % draw the plane if needed
    if config.collision_with_plane
      if config.collision_plane_z > min(verts(:,3))
%         config.collision_plane_z = min(verts(:,3)) - 20.0;
        config.collision_plane_z = min(verts(:,3));
      end
      hold on
      x_range = 1 * (max(verts(:,1))-min(verts(:,1)));
      y_range = 1 * (max(verts(:,2))-min(verts(:,2)));
      p = patch('XData', x_range*[-1 -1 1 1],'YData', y_range*[-1 1 1 -1], ...
                'ZData', config.collision_plane_z * [1 1 1 1], ...
                'FaceColor',[0.5 0.5 0.5], 'FaceAlpha', 0.5, 'FaceLighting', 'gouraud');
      axis equal
      h = plot3([], [], []);
      hold off
    end
    
    % collision with multiple sphere
    if config.collision_with_sphere
      [sphere_V, sphere_F] = readOBJ('./models/sphere.obj');
      if size(config.collision_sphere_c,1) < 1
        sphere_c = [mean(verts(:,1)) mean(verts(:,2)) max(verts(:,3)) + 0.5];
      else
        sphere_c = config.collision_sphere_c;
      end
      sphere_r = config.collision_sphere_r;
      sphere_V = sphere_r * sphere_V;
      sphere_vol = (4/3) * pi * sphere_r^3;
      sphere_v = repmat(config.collision_sphere_initial_veloity, size(sphere_c, 1), 1);
      sphere_rho = config.collision_sphere_rho;
      sphere_m = sphere_rho * sphere_vol;
      
      sphere_plt = cell(size(sphere_c, 1), 1);
      sphere_Vi = cell(size(sphere_c, 1), 1);
      sphere_face_normal = cell(size(sphere_c, 1), 1);
      hold on
      for si = 1:size(sphere_c, 1)
        sphere_Vi{si} = sphere_V + repmat(sphere_c(si,:), size(sphere_V,1), 1);
        sphere_plt{si} = tsurf(sphere_F, sphere_Vi{si});
        sphere_face_normal{si} = normals(sphere_Vi{si}, sphere_F);
      end
      hold off
    end
          
    collide_plane_vid = [];
    
    if config.save_obj
      mkdir(config.save_obj_path);
    end

    ii=1;
    for t=0:config.dt:40
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
        
        % Stiffness matrix (mex function)
        K = -vem3dmesh_neohookean_dq2(c, vol, params, dF_dc, w_I, k, n, ...
                                      size(x0_coms,2));
        K = L' * K * L;
%         K = J' * (P * (L' * K * L) * P') * J;
%           figure(3);spy(K);

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
%             dV_dF = SNH_dF(F, params(i,1), params(i,2));
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
          f_collision_other = zeros(3, 1);
          for i = 1:size(IF,1)
            fid = IF(i, 1); % colliding face id of the origin mesh
            c_fid = IF(i, 2); % colliding face id of the colliding mesh
            [f_collision, f_collision_other_i] = check_collision_between_faces(verts,faces,collide_verts,collide_faces,fid,c_fid,face_normal,collision_ratio,f_collision);
            f_collision_other = f_collision_other + f_collision_other_i;
          end
          
          if config.collision_with_other_sim
            % update the position of the other mesh
            collide_other_v = collide_other_v + config.dt * (f_collision_other + collide_other_m * config.gravity) /collide_other_m;
            collide_verts = collide_verts + config.dt * repmat(collide_other_v', size(collide_verts, 1), 1);
            % Always do a redraw
            t_collide.Vertices = collide_verts;
            drawnow
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
            [f_collision, ~] = check_collision_between_faces(verts,faces,verts,faces,fid,c_fid,face_normal,collision_ratio,f_collision);
            [f_collision, ~] = check_collision_between_faces(verts,faces,verts,faces,c_fid,fid,face_normal,collision_ratio,f_collision);
          end
        end
        % collision with the plane
        if config.collision_with_plane
          [collide_plane_vid, ~] = find(verts(:,3) < config.collision_plane_z);
          f_collision_plane_tmp = zeros(size(verts,1), 3);
          f_collision_plane_tmp(collide_plane_vid, 3) = 10 * collision_ratio * ...
                                                    (config.collision_plane_z*ones(size(collide_plane_vid,1), 1) - verts(collide_plane_vid,3));
          f_collision_plane = reshape(f_collision_plane_tmp', [], 1);                           
          f_collision = f_collision + f_collision_plane;
        end
        % collision with multiple sphere
        if config.collision_with_sphere
          for si = 1:size(sphere_c, 1)
            [IF] = intersect_other(verts,faces,sphere_Vi{si},sphere_F);
            f_collision_sphere_i = zeros(1, 3);
            for i = 1:size(IF,1)
              fid = IF(i, 1); % colliding face id of the origin mesh
              c_fid = IF(i, 2); % colliding face id of the colliding mesh              
              [f_collision, f_collision_sphere_tmp] = check_collision_between_faces(verts,faces,sphere_Vi{si},sphere_F,fid,c_fid,face_normal,collision_ratio,f_collision);
              f_collision_sphere_i = f_collision_sphere_i + f_collision_sphere_tmp'; 
            end
            % also check collision with other spheres
            dist_vec = repmat(sphere_c(si, :), size(sphere_v, 1), 1) - sphere_c;   
            dist_vec(si, :) = [];
            dist = sqrt(sum(dist_vec.^2, 2));
            [collide_idx, ~] = find(dist < sphere_r);
            if size(collide_idx, 1) > 0
              collide_idx
              dist_dir = dist_vec ./ repmat(dist, 1, 3);
              f_collision_other_sphere_i = config.collision_ratio * sum(repmat(2*sphere_r*ones(size(dist,1),1) - dist, 1, 3) .* dist_dir);
            else
              f_collision_other_sphere_i = zeros(1, 3);
            end
            % also check collision with the plane
            if config.collision_with_plane && sphere_c(si,3) < config.collision_plane_z + sphere_r
              f_collision_sphere_plane_i = collision_ratio * (config.collision_plane_z + sphere_r - sphere_c(si,3)) * [0 0 1];
            else
              f_collision_sphere_plane_i = zeros(1, 3);
            end
            % update velocity
            sphere_v(si, :) = sphere_v(si, :) + config.dt * ...
                              (f_collision_sphere_i + f_collision_other_sphere_i + f_collision_sphere_plane_i + sphere_m * [0 0 config.gravity] - sphere_v(si, :)) / sphere_m;
            % hacky air friction
            if size(IF,1) == 0 && sphere_v(si, 3) > 0
              if ii > 300
                sphere_v(si, 3) = 0.5 * sphere_v(si, 3);
              else
                sphere_v(si, 3) = 1.0 * sphere_v(si, 3);
              end
            elseif size(IF,1) > 0 && sphere_v(si, 3) < 0 && ii > 400
              sphere_v(si, 3) = 0;
            end
          end
          sphere_c = sphere_c + sphere_v * config.dt;
          
          % draw the spheres 
          hold on;
         	for si = 1:size(sphere_c, 1)
            sphere_Vi{si} = sphere_V + repmat(sphere_c(si,:), size(sphere_V,1), 1);
            set(sphere_plt{si}, 'Faces', sphere_F, 'Vertices', sphere_Vi{si});
            sphere_face_normal{si} = normals(sphere_Vi{si}, sphere_F);
          end
          hold off;
        end
        
        fprintf("f_internal = %d\n", sqrt(sum(f_internal.^2)));
        fprintf("f_gravity = %d\n", sqrt(sum(f_gravity.^2)));
        fprintf("f_error = %d\n", sqrt(sum(f_error.^2)));
        fprintf("f_collision = %d\n", sqrt(sum(f_collision.^2)));

        f_collision = hires_J' * f_collision;

        % Computing linearly-implicit velocity update
        lhs = J' * (P*(M - config.dt*config.dt*K)*P') * J;
        rhs = J' * (P*M*P'*J*qdot + f_internal + f_gravity + f_error + f_external) + f_collision;
        qdot = lhs \ rhs;
        
        if config.collision_with_plane && size(collide_plane_vid, 1) > 5 && ii >= 300
          qdot = reshape(qdot, 3, [])';
          qdot(:, [1 2]) = 0 * qdot(:, [1 2]);
          qdot = reshape(qdot', [], 1); 
        end
        
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
            obj_fn = config.save_obj_path + "part_" + int2str(ii) + ".obj";
            nurbs_write_obj(q,parts,obj_fn,ii);
            if config.collision_with_sphere
               for si = 1:size(sphere_c, 1)
                 obj_fn = config.save_obj_path + "sphere_" + int2str(si) + "_" + int2str(ii) + ".obj";
                 writeOBJ(obj_fn, sphere_Vi{si}, sphere_F);
               end
            end
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

function [f_collision, f_collision_other] = check_collision_between_faces(verts,faces,collide_verts,collide_faces,fid,c_fid,face_normal,collision_ratio,f_collision)
  % check vertex-face collisions in the colliding faces
  f_collision_other = zeros(3, 1);
  for j = 1:3  
    cv = collide_verts(collide_faces(c_fid,j),:);
    [dist, cp] = pointTriangleDistance(...
                            [verts(faces(fid,1),:); verts(faces(fid,2),:); verts(faces(fid,3),:)], ...
                            cv);
    if dot(cv-cp, face_normal(fid, :)) > 0
      f_collision(3*faces(fid,1)-2:3*faces(fid,1), 1) = f_collision(3*faces(fid,1)-2:3*faces(fid,1), 1) - 1/3 * collision_ratio * dist * face_normal(fid,:)';
      f_collision(3*faces(fid,2)-2:3*faces(fid,2), 1) = f_collision(3*faces(fid,2)-2:3*faces(fid,2), 1) - 1/3 * collision_ratio * dist * face_normal(fid,:)';
      f_collision(3*faces(fid,3)-2:3*faces(fid,3), 1) = f_collision(3*faces(fid,3)-2:3*faces(fid,3), 1) - 1/3 * collision_ratio * dist * face_normal(fid,:)';
      f_collision_other = f_collision_other + collision_ratio * dist * face_normal(fid,:)';
      disp('face collide!!');
    end
  end
end


