function z = ftrapz1d(x,y)
    x=x(:)';
    y=y(:);
    n = length(y);
    % Trapezoid sum computed with vector-matrix multiply.
    % NOTE: x and y need have different dimension below
    z = diff(x) * (y(1:n-1) + y(2:n))/2;
end % trapz1d
