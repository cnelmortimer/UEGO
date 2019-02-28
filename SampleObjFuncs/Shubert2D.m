function [fval] = Shubert2D(x)
    a=1.42513;
    b=0.80032;

    x1=x(1);
    x2=x(2);

    suma=0;
    for i=1:5
        suma=suma+i*cos((i+1)*x1+i);
    end

    sumb=0;
    for i=1:5
        sumb=sumb+i*cos((i+1)*x2+i);
    end
    
    fval=suma*sumb+1/2*power(x1+a,2)+1/2*power(x2+b,2);
end
% x(i) in [-10, 10]. Global minimum at x=(-1.42513,-0.80032) -> -186.7309
