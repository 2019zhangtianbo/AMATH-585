% AMATH 585 HW2 Problem 3
% Tianbo Zhang 1938501
% Create the plots for a and b

% Define the conditions: x, f(x)
x = linspace(-5, 5, 1000);
f_x = f(x);

% Plot the function f
plot(x, f_x, 'k-', 'LineWidth', 1.5, 'DisplayName', 'f(x)');
hold on;

% Calculate interpolating polynomials with formula xi = -5 + 10i/n
% Try n with 5, 10, 15, 20
for n = [5, 10, 15, 20]
    x_node = define_x(1, n);
    f_x_node = f(x_node);
    p_x = lagrange_interpolation(x_node, x, f_x_node);
    plot(x, p_x, 'LineWidth', 1.2, 'DisplayName', sprintf('n=%d', n));
    hold on;
end

% Add graph details
title('Interpolating Polynomials with Chebyshev nodes');
xlabel('x');
ylabel('f(x)');
legend show;
hold off;

g_x = f_x - p_x;
plot(x, g_x, 'LineWidth', 1.2, 'DisplayName', 'f(x)-p_{20}(x)');
title('Difference between f(x) and p_{20}(x)');
xlabel('x');
ylabel('f(x)');
legend show;
hold off;

function p_x = lagrange_interpolation(x_node, x, f_x_node)
    n = length(x_node);
    p_x = zeros(size(x));
    for i = 1:n
        l_x = ones(size(x));
        for j = [1:i-1, i+1:n]
            l_i = (x - x_node(j)) / (x_node(i) - x_node(j));
            l_x = l_x .* l_i;
        end
        size(l_x)
        p_x = p_x + f_x_node(i) * l_x;
    end
end

function x = define_x(option, n)
    if option == 0
        x = -5:10/n:5;
    else
        x = zeros(n);
        for i = 0 : n
            x(i+1) = 5 * cos(i * pi / n);
        end
    end
end

function f_x = f(x)
    f_x = 1 ./ (1 + x.^2);
end
