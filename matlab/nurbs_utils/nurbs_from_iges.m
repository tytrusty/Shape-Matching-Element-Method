function [nurbs_objects] = nurbs_from_iges(filename, res, rescale)
    if nargin < 3
        rescale=1;
    end
%     nurbs_object = {};
    data=iges2matlab(filename);

    %TODO change to parameters
    nrays = 8;
    samples_per_ray = 8;
    untrimmed_res = 10;

    % TODO: need to return weights of each sample on the surface.
    nurbs_objects = nurbs_sample(data, nrays, samples_per_ray, untrimmed_res);

    for i=1:numel(nurbs_objects)
        % Computed generalized coordinates
        srf = data{nurbs_objects{i}.surf_ptr};
        [p, J] = nurbs_coords(srf.nurbs, nurbs_objects{i}.UV);
        J = permute(J,[3 2 1]);

        nurbs_objects{i}.p = p;
        nurbs_objects{i}.pdot=zeros(size(p));
        nurbs_objects{i}.J = J;
        nurbs_objects{i}.x0 = squeeze(sum(J .* p,1));
        nurbs_objects{i}.srf = srf;

        T = nurbs_triangulate(nurbs_objects{i}, untrimmed_res);
        nurbs_objects{i}.T = T;
    end
end