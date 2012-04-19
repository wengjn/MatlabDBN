%plots demonstrating importance sampling
%http://www.cs.cmu.edu/~ggordon/MCMC/impsampplots.m

function impsampplots()
clear all
close all
clc

sampsize = 30; %sample size to use for examples

%some data -- a function that is big near 0
xs = -1:.01:1;
ys = foo(xs);
plot(xs, ys, 'LineWidth', 2);
set(gca, 'FontSize', 24);
fprintf('true integral %f\n', sum(ys)*2/length(ys))

%integrate by uniform sampling
figure
samp = rand(sampsize, 1)*2-1;
plot(xs, ys, 'LineWidth', 2);
hold on;
plot(samp, foo(samp), 'rx', 'MarkerSize', 15);
hold off;
set(gca, 'FontSize', 24);
fprintf('uniform sampling %f\n', 2*sum(foo(samp))/length(samp));

%integrate by importance sampling w/ normal proposal
figure
sig = .25;
samp = randn(sampsize, 1)*sig;
wts = exp(-samp.^2/(2*sig^2))/sqrt(2*pi*sig^2);
plot(xs, ys, 'LineWidth', 2);
hold on;
plot(samp, foo(samp), 'rx', 'MarkerSize', 15);
hold off;
set(gca, 'FontSize', 24);
fprintf('importance sampling %f\n', sum(foo(samp)./wts)/length(samp));

%plot that demonstrates reason for bias of parallel IS
figure
plot(.3:.1:3, 1./(.3:.1:3), '-', [.5 2], [2 .5], 'x', [1.25 1.25], ...
    [.8 1.25], 'o', 'LineWidth', 3, 'MarkerSize', 15);
hold on;
plot([.5 2], [2 .5], '--');
hold off;
axis([0 3 0 3]);
set(gca, 'FontSize', 24);
xlabel('mean(weights)');
ylabel('1 / mean(weights)');
text(1.3, 2.5, 'E(mean(weights))', 'FontSize', 20);
line([1.25 1.25], [0 3], 'LineStyle', ':');
line(.5, 0, 'Marker', 'x', 'LineWidth', 3, 'MarkerSize', 15, 'Color', [0 .5 0])
line(2, 0, 'Marker', 'x', 'LineWidth', 3, 'MarkerSize', 15, 'Color', [0 .5 0])
line(1.25, 0, 'Marker', 'o', 'LineWidth', 3, 'MarkerSize', 15, 'Color', [1 0 0]);

%parallel IS example: we observe that we are facing a known 
%landmark and 0.9m away from it, but don't know which direction
%we are facing, so our posterior is shaped like part of a ring

n = 500; %sample size
ax = [-2 2 -2 2];  %plot region
xysig = 1;  %prior stdev of x, y
thsig = 10; %prior stdev of theta
osig = .3333  %observation stdev
obs = 0.9;   %observed distance to landmark
ox = 1;   %landmark x
oy = .3;  %landmark y

%pick samples from a normal proposal distribution
xs = xysig * randn(n, 1);
ys = thsig * randn(n, 1);
ths = thsig * randn(n, 1);

%plot the proposal samples and the ring around the landmark
figure
d = 0:2*pi/100:2*pi;
plot(xs, ys, 'b.', ox, oy, 'rx', ox+obs*cos(d), oy+obs*sin(d), 'r-', ...
    'LineWidth', 2);
axis equal;
axis(ax);

%compute importance weights
ws = sqrt((xs+obs*cos(ths)-ox).^2 + (ys+obs*sin(ths)-oy).^2);
ws = exp(-ws.^2/(2*osig*osig));

%plot importance weights
maxw = max(ws);
ws = ws./maxw;
plot(ox, oy, 'rx', ox+obs*cos(d), oy+obs*sin(d), 'r-', 'LineWidth', 2);
axis equal;
axis(ax);
for i = 1:n
    line(xs(i), ys(i), 'Marker', 'o', 'MarkerSize', 1+10*sqrt(ws(i)))
    if(ws(i)>0.5)
        line([xs(i), xs(i)+.1*cos(ths(i))], ...
            [ys(i), ys(i)+.1*sin(ths(i))]);
    end
end

%normalize theta to -pi, pi
while (sum(ths > pi) > 0 )
    ths(ths > pi) = ths(ths > pi) -2 * pi;
end

while (sum(ths < -pi) > 0)
    ths(ths < - pi) = ths(ths < -pi) + 2 * pi;
end
%display posterior mean
fprintf('Expected x y theta pos %f %f %f\n', (ws./sum(ws))' * [xs ys ths]);

%the function that's integrated in the importance sampling examples
function ys = foo(xs)
ys = sin(xs*8).^2.*(xs).^-2;
ys(xs==0) = 64;
return 