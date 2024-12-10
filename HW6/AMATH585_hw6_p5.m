% AMATH 585 HW6 Problem 5
% Tianbo Zhang 1938501
% Solve the IVP problem with Forward Euler, Classiscal Fourth-order
% Runge-kutta method

% Set up base conditions
T = 8;
y0 = 1;
N_values = [25, 50, 100, 200, 400, 800, 1600];
error_euler = zeros(1, length(N_values));
error_runge = zeros(1, length(N_values));

for i = 1:length(N_values)
    N = N_values(i);
    h = T/N;
    t = 0:h:T;

    %Forward Euler Method
    y_euler = zeros(1, N+1);
    y_euler(1) = y0;
    for j = 1:N
        y_euler(j+1) = y_euler(j) + h * dy_dt(t(j), y_euler(j));
    end

    % Classical Fourth-order Runge-kutta method
    y_runge = zeros(1, N+1);
    y_runge(1) = y0;
    for j = 1:N
        q1 = dy_dt(t(j), y_runge(j));
        q2 = dy_dt(t(j)+h/2, y_runge(j) + h*q1/2);
        q3 = dy_dt(t(j)+h/2, y_runge(j) + h*q2/2);
        q4 = dy_dt(t(j)+h, y_runge(j) + h*q3);
        y_runge(j+1) = y_runge(j) + (h/6)*(q1 + 2 * q2 + 2 * q3 + q4);
    end

    % Compute Errors
    y_true = y(t);
    error_euler(i) = max(abs(y_true - y_euler));
    error_runge(i) = max(abs(y_true - y_runge));

    % Plot for N+25
    if N == 25
        figure;
        plot(t, y_true, 'k-');
        hold on
        plot(t, y_euler, 'ro');
        hold on
        plot(t, y_runge, 'bo');
        hold off
        legend('True Solution', 'Euler Approximation', 'Runge Approximation');
        title('Solution for N=25');
        xlabel('Time t');
        ylabel('Solution y(t)');
    end
end

% Plot errors
h_values = T./N_values;
figure;
loglog(h_values, error_euler, 'r-o', h_values, error_runge, 'b-o');
legend('Error Euler', 'Error Runge');
title('Error vs. Step Size');
xlabel('Step Size h');
ylabel('Error');
grid on;

% Calculate y'
function dy = dy_dt(t, y)
    dy = y^2 - sin(t) - cos(t)^2;
end

% Calculate true solution
function y_true = y(t)
    y_true = cos(t);
end