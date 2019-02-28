function [fval]=Schwefel2D(x)
    fval=sum(-x.*sin(sqrt(abs(x))));
end
% x(i) in [-500, 500]. Global minimum at x=(420.9687, 420.9687) -> -837.9658
