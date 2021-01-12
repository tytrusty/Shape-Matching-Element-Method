function k = basis_size(d, order)
    k=0;
    % Iteration starts at 1 (rather than 0) because our monomial basis
    % excludes the constant term (we globally solve for 1 constant term)
    for i=1:order
        k=k+nchoosek(d+i-1,i);
    end
end