function dF_dc = vem_dF_dc(dM_dX, W, w, E, d, k)
    m = size(w,1);
    n = size(w,2);
    
    dF_dc = zeros(m, d*d, d*(k*n + 1));
    for i = 1:m
        dMi_dX = squeeze(dM_dX(i,:,:));
        Wi = squeeze(W(i,:,:));
        dF_dc(i,:,:) = dMi_dX * Wi;
    end
end