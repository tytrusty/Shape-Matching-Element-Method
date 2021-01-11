function [normals, areas] = nurbs_normals(srf, UV, p)
    [f, f_grad, g, g_grad] = nurbs_grad(srf, UV);

    n = size(UV,2);
    areas = zeros(n, 1);
    normals = zeros(n, 3);
    
    % What was i thinking when i wrote this
    for i=1:n
        fij = squeeze(f(i,:,:)) * p;
        fij_ugrad = squeeze(f_grad(i,1,:,:)) * p;
        fij_vgrad = squeeze(f_grad(i,2,:,:)) * p;
        g2 = g(i)*g(i);
        u_jac = (g(i) * fij_ugrad - fij * g_grad(i,1)) / g2;
        v_jac = (g(i) * fij_vgrad - fij * g_grad(i,2)) / g2;
        jac = [u_jac v_jac];
        normals(i,:) = cross(u_jac, v_jac);
        areas(i) = sqrt(det(jac'*jac));
    end
    normals = normals ./ vecnorm(normals,2,2);

    % TODO use this version of computing nurbs gradients. It's much faster
    %      than my crummy version.
    %     dcrv = nrbderiv(srf);
    %     [pnts,jac] = nrbdeval(srf, dcrv, UV);
    %     C = cross(jac{1}, jac{2});
    %     lengths = vecnorm(C,2,1)';
    %     normals = C ./ lengths;
end