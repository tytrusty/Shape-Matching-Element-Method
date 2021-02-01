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

%quadratic monomial basis without constant term
monomial = @(X) [X(1).*X(1) X(2).*X(2) X(1).*X(2) X(1) X(2) 0 0 0 0 0; ...  
                 0 0 0 0 0 X(1).*X(1) X(2).*X(2) X(1).*X(2) X(1) X(2)];

             
q = igl2bart(V); 

for jj = 1:50
    
q(4) = q(4) + 0.1;
    

%build projection hiearchachly by fitting low order terms first

%(1) shape match with constant term
%constant term normally would be ome point inside shape but I only have a
%line so I'll just set it to be the origin for now
a = [0;0];


%build monomial matrix
A = [monomial(V(1,:)- a'); monomial(V(2,:) - a'); monomial(V(3,:)- a')];
b = q;

%this is gross
A_hires = [];
for ii=1:size(V_hires,1) 
    A_hires = [A_hires; monomial(V_hires(ii,:) - a')];
end

%(2) shape match linear term with left overs
c_linear = (A(:, [4 5 9 10])'*A(:, [4 5 9 10]))\(A(:, [4 5 9 10])'*(b-repmat(a,3,1)));

%(3) shape match quadratic term with left overs
c_quadratic = (A(:, [1 2 3 6 7 8])'*A(:, [1 2 3 6 7 8]))\(A(:, [1 2 3 6 7 8])'*(b-repmat(a,3,1) - A(:, [4 5 9 10])*c_linear));

q_out = A_hires(:, [1 2 3 6 7 8])*c_quadratic + A_hires(:, [4 5 9 10])*c_linear + repmat(a,size(V_hires,1),1);
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