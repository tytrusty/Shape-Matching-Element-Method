function [normals, areas] = nurbs_normals(nrb, u, v, p, f, f_grad, g, g_grad)
    if nargin < 8
        [f, f_grad, g, g_grad] = nurbs_grad(nrb, u, v);
    end
    
    areas = zeros([size(u,1) size(v,1)]);
    normals = zeros([size(u,1) size(v,1) 3]);
    for i=1:size(u,1)
        for j=1:size(v,1)
            fij = squeeze(f(i,j,:,:)) * p;
            fij_ugrad = squeeze(f_grad(i,j,1,:,:)) * p;
            fij_vgrad = squeeze(f_grad(i,j,2,:,:)) * p;
            g2 = g(i,j)*g(i,j);
            u_jac = (g(i,j) * fij_ugrad - fij * g_grad(i,j,1)) / g2;
            v_jac = (g(i,j) * fij_vgrad - fij * g_grad(i,j,2)) / g2;
            jac = [u_jac v_jac];
            normals(i,j,:) = cross(u_jac, v_jac);
            
            
            areas(i,j) = sqrt(det(jac'*jac));
       end
    end
    % TODO: weighted normals & areas are ridiculous...
    % look into more.
    normals = normals ./ vecnorm(normals,2,3);

end