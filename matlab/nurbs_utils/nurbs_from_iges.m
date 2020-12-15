function [nurbs_object] = nurbs_from_iges(filename, res, rescale)
    if nargin < 3
        rescale=1;
    end
    nurbs_object = {};
    data=iges2matlab(filename);
%     [data,EntityType,numEntityType,unknownEntityType,numunknownEntityType]=iges2matlab(filename)
    ii=1;
    for i = 1:numel(data)
        if data{i}.type == 128
            nrb = data{i}.nurbs;
            
            ratio = data{i}.ratio;
            
            if ~rescale
               ratio = [1 1]; 
            end
            
            if numel(res) > 1
                res_ii = res(ii);
            else
                res_ii = res;
            end
            
            subd = max(round(res_ii .*ratio ./ max(ratio)), 5);
            u_range = data{i}.u;
            v_range = data{i}.v ;

            % gauss-legendre positions & weights
            %[u,w_u]=lgwt(res,u_range(1),u_range(2));
            %[v,w_v]=lgwt(res,v_range(1),v_range(2));

            % trapezoidal quadrature
            u = linspace(u_range(1),u_range(2),subd(1))';
            v = linspace(v_range(1),v_range(2),subd(2))';
            w_u = ones(subd(1),1) * (u_range(2)-u_range(1)) / subd(1) / 2;
            w_v = ones(subd(2),1) * (v_range(2)-v_range(1)) / subd(2) / 2;
            % comment out on periodic surface
            w_u(2:end-1) = w_u(2:end-1) * 2;
            w_v(2:end-1) = w_v(2:end-1) * 2;

            % Computed generalized coordinates
            [p, J] = nurbs_coords(data{i}.nurbs, u, v);
            J_flat=reshape(J,[],3,size(J,4));
            J_flat=permute(J_flat,[3 2 1]);

            surf.u_range = u_range;
            surf.v_range = v_range;
            surf.u = u;
            surf.v = v;
            surf.w_u = w_u;
            surf.w_v = w_v;
            surf.p = p;
            surf.pdot=zeros(size(p));
            surf.J = J;
            surf.J_flat = J_flat;
            surf.x0 = squeeze(sum(J_flat .* p,1));
            surf.res = res;
            surf.subd = subd;
            surf.nurbs = data{i}.nurbs;
            nurbs_object{ii} = surf;
            ii = ii+1;
        end
    end
end