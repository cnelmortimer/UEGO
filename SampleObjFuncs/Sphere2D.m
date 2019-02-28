function [feval] = Sphere2D(x)
    feval = x(1)^2 + x(2)^2;
end
% See: https://en.wikipedia.org/wiki/Test_functions_for_optimization
% x(i) in [-inf, inf]. Global minimum at x=(0, 0) -> 0
