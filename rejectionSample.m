%rejection sampling algorithm to compute the integration of y = xexp(x)dx
%[0 1]

%f(x) is the PDF of x, h(x) is the proposal PDF, we choose as uniform
%distribution, since f(x) is less than 1.6, we choose c = 1.6
%we have f(x) <= c*h(x)

% N = 1000;
% for i = 1:N
%     u = 2;
%     alpha = 1;
%     while (u > alpha)
%         xminus = rand;
%         u = rand;
%         alpha = exp(xminus)/(exp(1)-1)/(1.6*u);
%     end
%     x(i) = xminus;
% end
% 
% y = (exp(1) - 1)*sum(x) / N;


N = 1000;
for i = 1:N
    u = 2;
    alpha = 1;
    while (u > alpha)
        xminus = rand;
        u = rand;
        alpha = exp(xminus) / (1.6 * u);
    end
    x(i) = xminus;
end

y = (exp(1) -1)*sum(x) / N;
