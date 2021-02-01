function vem_2d_test_matching
    % Simulation parameters
    order = 1;          % (1 or 2) linear or quadratic deformation
	d = 2;              % dimension (2 or 3)
    save_output = 0;    % (0 or 1) whether to output images of simulation
    
    % The number of elements in the monomial basis.
    k = basis_size(1, order);
    
    % Shape matching modes:
%     mode = 'global';               % global solve with regular inverse
%     mode = 'global_pinv';           % global solve with pseudoinverse
%     mode = 'global_svd_truncated';  % global solve with truncated svd inv
%     mode = 'local';
    mode = 'local_pinv';
%     mode = 'local_svd_truncated';
    
    offset = 0;
    shape_func = @(a,n) ((a - offset).^n);
    
    range = [0 1];
    samples = [20 30];
    n = 1;
    n_end = 5;
    x = linspace(range(1), range(2), samples(1));
    y = shape_func(x,n);
    
    y_com = mean(y);
    L = compute_shape_matrices(y, y_com, {1:samples(1)}, order, mode);
    
    c = L * y(:);
    
    x0 = linspace(range(1), range(2), samples(2));
    x0_com = mean(x0);
    
    M = monomial_basis(x0, x0_com, order);
    
    y0 = c(1:k)' * M + c(end);
    
    clf;
    plot(x,y);
    hold on;
    plot(x0,y0,'.','MarkerSize',20);
    axis equal;
    
    deformations=100;
    n_vals = linspace(n,n_end, deformations);
    
    
    for i=1:deformations
        clf;
        y_new = shape_func(x, n_vals(i));
        c = L * y_new(:);
        y_pred = c(1:k)' * M + c(end);
        plot(x,y_new);
        hold on;
        plot(x0,y_pred,'.','MarkerSize',20);
        axis equal;
        legend('Deformed Positions', 'Predicted positions','Location','northwest');
        drawnow;
    end
       
end