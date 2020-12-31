function [nurbs_objects] = nurbs_from_iges(filename, res, rescale)
    if nargin < 3
        rescale=1;
    end
%     nurbs_object = {};
    data=iges2matlab(filename);

    % TODO: need to return weights of each sample on the surface.
    enable_trimming = 0;
    use_triangle = 0;
    nurbs_objects = nurbs_sample(data, enable_trimming);
    disp('Done Generating samples');

    for i=1:numel(nurbs_objects)
        % Computed generalized coordinates
        srf = data{nurbs_objects{i}.surf_ptr};
        nurbs_objects{i}.srf = srf;
                
        [p, J] = nurbs_coords(srf.nurbs, nurbs_objects{i}.UV);
        J = permute(J,[3 2 1]);

        nurbs_objects{i}.p = p;
        nurbs_objects{i}.pdot=zeros(size(p));
        nurbs_objects{i}.J = J;
        nurbs_objects{i}.x0 = squeeze(sum(J .* p,1));

        [T,TUV] = nurbs_triangulate(nurbs_objects{i}, use_triangle);
        nurbs_objects{i}.T = T;
        
        % Constructing jacobian for Trimesh UV coordinates
        [~, T_J] = nurbs_coords(srf.nurbs, TUV);
        T_J = permute(T_J, [3 2 1]);
        nurbs_objects{i}.T_J = T_J;
        nurbs_objects{i}.Vertices = squeeze(sum(T_J .* p,1));
    end
end