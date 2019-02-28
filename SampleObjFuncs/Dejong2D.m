function fval = Dejong2D(x)
    fval = (x(1))^2 + (x(1)+x(2))^2;
end
% See: https://en.wikipedia.org/wiki/Test_functions_for_optimization (Schwefel's function)
% x(i) in [-5.12, 5.12]. Global minimum at x=(0, 0) -> 0
