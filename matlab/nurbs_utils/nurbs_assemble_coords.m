function [J,q,E,x] = nurbs_assemble_coords(parts)

    x = zeros(3,0);
    q_size = 0;
    J_size = 0;

    % build global position vectors
    for i=1:numel(parts)
        idx1=q_size+1;
    	q_size = q_size + size(parts{i}.p,1);
        J_size = J_size + size(parts{i}.J,2) * size(parts{i}.J,3);
        idx2=q_size;
        
        % indices into global configuration vector
        parts{i}.idx1=idx1;
        parts{i}.idx2=idx2;
    end
    
    q = zeros(q_size,1);
    J = zeros(J_size,q_size);
    J_idx = [0 0];
    for i=1:numel(parts)
        q(parts{i}.idx1:parts{i}.idx2,1) = parts{i}.p;
        Ji = parts{i}.J(:,:)';
        cond(Ji'*Ji);
        % Block indices
        I1=J_idx(1)+1:J_idx(1)+size(Ji,1);
        I2=J_idx(2)+1:J_idx(2)+size(Ji,2);
        
        % Subsitute block into global NURBs jacobian (J) matrix
        J(I1,I2) = Ji;
        J_idx = J_idx + size(Ji);
    end
    J = sparse(J);
    E=cell(numel(parts),1);
    
    for i=1:numel(parts)
        idx1=size(x,2)+1;
        x = [x parts{i}.x0];
        idx2=size(x,2);
        
        % indices into global configuration vector
        E{i}=idx1:idx2;
    end

end