function vem_fitting_test2
%polynomial fitting test
%define quadratic projection for 3 point line
%idea: shape match should minimize deformation
%conjecture: corresponsds to keeping higher order terms as close to zero
%as possible
%requirement: surface patch must at least be able to uniquely determine 
% the linear term (deformation gradient)
% this requires explcitly setting the constant term as we've been doing

V = [-1 1; 0 1; 1 1];

V_hires = [linspace(-1,1, 100)' linspace(1,1, 100)'];      
q = igl2bart(V); 

COM = [0 0.9];

L = compute_shape_matrices(V', COM', {[1 2 3]}, 2, 'hierarchical');
Y = monomial_basis_matrix(V_hires',  COM', 2, 5);

for jj = 1:50

q(3) = q(3) + 0.01;
q(4) = q(4) + 0.1;
    
c = L*(q);

q_outx = squeeze(Y(:,1,:)) * c(1:end-2) + c(end-1);
q_outy = squeeze(Y(:,2,:)) * c(1:end-2) + c(end);
q_out = [q_outx'; q_outy'];
q_out = q_out(:);

plot(q_out(1:2:end), q_out(2:2:end));
axis([-1 1 0 5]);
hold on
    plot(q_out(1:2:end), q_out(2:2:end), 'b*');
    
    plot(q(1:2:end), q(2:2:end), 'g');
    plot(q(1:2:end), q(2:2:end), 'g*');
hold off

pause(0.1)
drawnow

end
end