function [nurbs_objects] = nurbs_from_iges(filename, sample_density)
    if nargin < 2
        sample_density=2;
    end
%     nurbs_object = {};
    data=iges2matlab(filename);

    % TODO: need to return weights of each sample on the surface.
    enable_trimming = 1;
    use_triangle = 1;
    generate_hires = 1;
    trimmed_hires_res = 8;
    nurbs_objects = nurbs_sample(data, enable_trimming, sample_density);
    disp('Done Generating samples');

    for i=1:numel(nurbs_objects)
        % Computed generalized coordinates
        srf = data{nurbs_objects{i}.surf_ptr};
        nurbs_objects{i}.srf = srf;
                
        [p, J] = nurbs_coords(srf.nurbs, nurbs_objects{i}.UV);
        J = permute(J,[3 2 1]);

        Ji = J(:,:)';
        JJ = Ji'*Ji;
        cond(JJ)

        nurbs_objects{i}.p = p;
        nurbs_objects{i}.pdot=zeros(size(p));
        nurbs_objects{i}.J = J;
        nurbs_objects{i}.x0 = squeeze(sum(J .* p,1));

        [T,~] = nurbs_triangulate(nurbs_objects{i});
        nurbs_objects{i}.T = T;
        
        if generate_hires
            [T,TUV] = nurbs_triangulate(nurbs_objects{i}, ...
                                        use_triangle, trimmed_hires_res);
            [~, J] = nurbs_coords(srf.nurbs, TUV);
            J = permute(J,[3 2 1]);
            nurbs_objects{i}.hires_T = T;
            nurbs_objects{i}.hires_J = J;
            nurbs_objects{i}.hires_x0 = squeeze(sum(J .* p,1));
            nurbs_objects{i}.hires_UV = TUV;
        end
    end
end