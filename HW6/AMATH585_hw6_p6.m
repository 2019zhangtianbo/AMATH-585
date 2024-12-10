% AMATH 585 HW6 Problem 6
% Tianbo Zhang 1938501
% Use ode45 to solve the Lotka-Volterra predator-prey equations

% Set up initial equation
R0 = 20;
F0 = 20;
y0 = [R0; F0];

% Define time span
t_span = [0, 50]

[T, Y] = ode45(@lotka_volterra, t_span, y0);

% Plot R(t) and F(t) vs. t
figure;
plot(T, Y(:,1), 'r-', T, Y(:,2), 'b-');
xlabel('Time');
ylabel('Population');
legend('Prey (R)', 'Predator (F)');
title('Prey and Predator populations over time');

% Plot F vs. R
figure;
plot(Y(:,1), Y(:,2), 'm-');
xlabel('Prey Population (R)');
ylabel('Predator Population (F)');title('Predator vs. Prey Population Dynamics');

function dydt = lotka_volterra(t, y)
    dydt =[(1 - 0.02 * y(2)) * y(1); (-1 + 0.03 * y(1)) * y(2)];  
end