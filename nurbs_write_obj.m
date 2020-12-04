function nurbs_write_obj(q, parts, res, path)
    faces=[];
    verts=[];
    for i=1:numel(parts)
        obj_u = linspace(parts{i}.u_range(1),parts{i}.u_range(2),res)';
        obj_v = linspace(parts{i}.v_range(1),parts{i}.v_range(2),res)';
        [~, obj_J] = nurbs_coords(parts{i}.nurbs, obj_u, obj_v);
        obj_J = reshape(obj_J,[],3,size(obj_J,4));
        obj_J = permute(obj_J,[3 2 1]);
        obj_q = q(parts{i}.idx1:parts{i}.idx2);
        obj_x = squeeze(sum(obj_J .* obj_q,1));
        obj_x = reshape(obj_x, 3, res, res);
        
        % Create surface mesh from point samples
        plt = surf( squeeze(obj_x(1,:,:)), ...
                    squeeze(obj_x(2,:,:)), ...
                    squeeze(obj_x(3,:,:)));
                
        % Convert surface to patch to give us a vertex & face set.
        fvc = surf2patch(plt);
        delete(plt);
        fvc.faces = fvc.faces + size(verts,1);
        faces=[faces; fvc.faces]; % sue me
        verts=[verts; fvc.vertices];
    end
    
    % you better have gptoolbox, son
    writeOBJ(path, verts, faces);
end

