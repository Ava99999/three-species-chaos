function X = superpredators(t,x,f1,f2,d1,d2)
    X    = zeros(3,1);
    X(1) = x(1)*(1-x(1)) - f1(x(1))*x(2);
    X(2) = f1(x(1))*x(2) - f2(x(2))*x(3) - d1*x(2);
    X(3) = f2(x(2))*x(3) - d2*x(3);
end