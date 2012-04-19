%importance sampling to compute integration of y = xexp(x)dx

%f(x) is the PDF of x, h(x) is the proposal PDF
%wherever f(x)>0, h(x)>0

N = 1000;
for i = 1:N
    while 1 > 0
        r = 1 + 0.25 * randn;
        if r>=0 && r<=1
            break;
        end
    end
    
    %r = 1 + 0.25 * randn;
    x(i) = r;
    w(i) = sqrt(2*pi) * exp(x(i) + 2*(x(i) - 1)^2) / (4*(exp(1) - 1));
end

y = (exp(1) - 1) * sum(x .* w) / N;