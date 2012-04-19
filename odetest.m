function xprime = odetest(t, x)
%x(1) = x, x(2) = y, x(3) = z

xprime(1) = -x(2)*exp(-t/5)+x(3)*exp(-t/5)+1;
xprime(2) = x(3);
xprime(3) = -2*sin(t);

xprime = xprime(:);