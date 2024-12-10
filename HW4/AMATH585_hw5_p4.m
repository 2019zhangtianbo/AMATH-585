

T = 2 * pi;
alpha = 0.7;
beta = 0.7;
n = 10;
h = T / n;
t = linspace(0, T, n+1);
theta = linspace(0, 2*pi, n+1);
theta = theta(2:end-1)
% theta = 0.7 * cos(t(2:end-1)) + 0.5 * sin(t(2:end-1));
size(theta)
plot(t(2:end-1), theta, 'DisplayName', sprintf('Iteration 0'));
hold on;
for iteration = 1:10
    F = getF(n, h, theta, alpha, beta);
    J = getJacobian(n, h, theta);
    delta = J\(-F);
    theta = theta + delta';
    plot(t(2:end-1), theta, 'DisplayName', sprintf('Iteration %d', iteration));
    hold on;
end
hold off;
xlabel('t');
ylabel('\theta(t)');
title('Approximate Solution with Linear Interpolation');
legend show;

function J = getJacobian(n, h, theta)
    J = zeros(n-1, n-1);
    J(1, 1:2) = [-2/h^2 + cos(theta(1)), 1/h^2];
    for i = 2:n-2
        J(i, i-1:i+1) = [1/h^2, -2/h^2 + cos(theta(i)), 1/h^2];
    end
    J(n-1, n-2:n-1) = [1/h^2, -2/h^2 + cos(theta(n-1))];
end

function F = getF(n, h, theta, alpha, beta)
    F = zeros(n-1, 1);
    F(1) = (alpha - 2*theta(1) + theta(2)) / h^2 + sin(theta(1));
    for i = 2:n-2
        F(i) = (theta(i-1) - 2*theta(i) + theta(i+1)) / h^2 + sin(theta(i));
    end
    F(n-1) = (theta(n-2) - 2*theta(n-1) + beta) / h^2 + sin(theta(n-1));
end