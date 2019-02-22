function f_x = Rastrigin2D(x)
	f_x = 20.0 + x(1)^2 + x(2)^2 - 10.0*(cos(2*pi*x(1)) + cos(2*pi*x(2)));
end
% x(i) in [-5.12, 5.12]. Global minimum at x=(0, 0) -> 0
% See: https://es.mathworks.com/help/gads/example-rastrigins-function.html
% See: https://en.wikipedia.org/wiki/Rastrigin_function
