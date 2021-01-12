function dF_dc = vem_dF_dc(dM_dX, W)
    m = size(dM_dX,1);
    dF_dc = cell(m,1);
    
    for i = 1:m
        dMi_dX = squeeze(dM_dX(i,:,:));
        dF_dc{i} = dMi_dX * W{i};
    end
end