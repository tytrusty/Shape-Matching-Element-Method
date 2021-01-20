function nurbs_write_iges(q, parts, path)
    
    nurbs_objects=cell(numel(parts),1);
    q_idx = 0;
    for i=1:numel(parts)
        srf = parts{i}.srf;
        
        qi = q(q_idx+1:q_idx+numel(parts{i}.p));
        qi = reshape(qi, [3 size(srf.nurbs.coefs,3) size(srf.nurbs.coefs,2)]);
        qi = permute(qi, [1 3 2]);
    
        
        srf.q = qi(:);
        nurbs_objects{i} = srf;
        q_idx = q_idx + numel(parts{i}.p);
    end
    
   igesout(nurbs_objects, path)
end

