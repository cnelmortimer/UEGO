function [fval] = Rosenbrock2D(x)
    x1=x(1);
    x2=x(2);
    fval=0.5*(x1.^2-x2).^2+(x1-1).^2;
end
% x(i) in [-3, 3]. Global minimum at x=(1, 1) -> 0
% See: https://es.wikipedia.org/wiki/Funci?n_de_Rosenbrock
