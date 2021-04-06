function [J,hires_J,q,E,x,S] = nurbs_assemble_coords(parts)

    x = zeros(3,0);
    q_size = 0;
    J_size = 0;
    hires_J_size = 0;
    S=[];

    % build global position vectors
    for i=1:numel(parts)
        idx1=q_size+1;
      	q_size = q_size + size(parts{i}.p,1);
        J_size = J_size + size(parts{i}.J,2) * size(parts{i}.J,3);
        hires_J_size = hires_J_size + size(parts{i}.hires_J,2) * size(parts{i}.hires_J,3);
        idx2=q_size;
        
        % need to repeat each weight 3x t account for d=3
        Si = [parts{i}.W parts{i}.W parts{i}.W]';
        S = [S; Si(:)];

        % indices into global configuration vector
        parts{i}.idx1=idx1;
        parts{i}.idx2=idx2;
    end
    
    % create sparse diagonal weight matrix
    S=spdiags(S,0,numel(S),numel(S));

    q = zeros(q_size,1);
    J = zeros(J_size,q_size);
    hires_J = zeros(hires_J_size,q_size);
    J_idx = [0 0];
    hires_J_idx = [0 0];
    for i=1:numel(parts)
        q(parts{i}.idx1:parts{i}.idx2,1) = parts{i}.p;
        Ji = parts{i}.J(:,:)';
        hires_Ji = parts{i}.hires_J(:,:)';

        % Block indices
        I1=J_idx(1)+1:J_idx(1)+size(Ji,1);
        I2=J_idx(2)+1:J_idx(2)+size(Ji,2);
        
        hires_I1=hires_J_idx(1)+1:hires_J_idx(1)+size(hires_Ji,1);
        hires_I2=hires_J_idx(2)+1:hires_J_idx(2)+size(hires_Ji,2);
        
        % Subsitute block into global NURBs jacobian (J) matrix
        J(I1,I2) = Ji;
        J_idx = J_idx + size(Ji);
        
        hires_J(hires_I1,hires_I2) = hires_Ji;
        hires_J_idx = hires_J_idx + size(hires_Ji);
    end
    J = sparse(J);
    hires_J = sparse(hires_J);
    E=cell(numel(parts),1);
    
    for i=1:numel(parts)
        idx1=size(x,2)+1;
        x = [x parts{i}.x0];
        idx2=size(x,2);
        
        % indices into global configuration vector
        E{i}=[idx1:idx2]';
    end

end