function nurbs_write_obj(q, parts, path, ii)
    faces=[];
    verts=[];

    q_idx = 0;
    for i=1:numel(parts)
        
        qi = q(q_idx+1:q_idx+numel(parts{i}.p));
        xi = squeeze(sum(parts{i}.hires_J .* qi,1));

        faces_i = parts{i}.hires_T + size(verts,1);
        faces=[faces; faces_i];
        verts=[verts; xi'];
        
        % obj_fn = "output/obj/part_" + i + "_" + int2str(ii) + ".obj";
        % writeOBJ(obj_fn, fvc.vertices, fvc.faces);
        q_idx = q_idx + numel(parts{i}.p);
    end
    
    % you better have gptoolbox, son
    writeOBJ(path, verts, faces);
end

