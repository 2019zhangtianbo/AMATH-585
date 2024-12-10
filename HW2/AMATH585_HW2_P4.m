% AMATH 585 HW2 Problem 4 (Textbook p125, problem27)
% Tianbo Zhang 1938501
% Create the plots for a and b

f1 = figure;
f2 = figure;
f3 = figure;
for n = [5, 10, 20]
    x_a = define_xa(n);
    x_b = define_xb(n);
    x_c = define_xc(n);
    grid_a = found_grid(x_a);
    grid_b = found_grid(x_b);
    grid_c = found_grid(x_c);
    lamb_a = lagrange_interpolation(x_a, grid_a);
    lamb_b = lagrange_interpolation(x_b, grid_b);
    lamb_c = lagrange_interpolation(x_c, grid_c);
    
    delta_a = max(lamb_a);
    delta_b = max(lamb_b);
    delta_c = max(lamb_c);
    fprintf('(a)The Lebesgue constant for n = %d is %d \n', n, delta_a);
    fprintf('(b)The Lebesgue constant for n = %d is %d \n', n, delta_b);
    fprintf('(c)The Lebesgue constant for n = %d is %d \n', n, delta_c);

    figure(f1);
    plot(grid_a, log10(lamb_a), 'LineWidth', 1.2, 'DisplayName', sprintf('n=%d', n));
    hold on;

    figure(f2);
    plot(grid_b, lamb_b, 'LineWidth', 1.2, 'DisplayName', sprintf('n=%d', n));
    hold on;

    figure(f3);
    plot(grid_c, lamb_c, 'LineWidth', 1.2, 'DisplayName', sprintf('n=%d', n));
    hold on;
end
figure(f1);
title('Lebsgue Function with equidistant nodes');
xlabel('x');
ylabel('f(x)');
legend show;
hold off;

figure(f2);
title('Lebsgue Function with Chebyshev nodes');
xlabel('x');
ylabel('f(x)');
legend show;
hold off;

figure(f3);
title('Lebsgue Function with Chebyshev nodes HW3 version');
xlabel('x');
ylabel('f(x)');
legend show;
hold off;

function grid_x = found_grid(x_n)
    grid_x = [x_n(1)];
    for i = 1:length(x_n)-1
        sub_int = linspace(x_n(i), x_n(i+1), 21);
        grid_x = [grid_x, sub_int(2:end)];
    end
end

function lambda = lagrange_interpolation(x_node, x)
    n = length(x_node);
    lambda = zeros(size(x));
    for i = 1:n
        l_x = ones(size(x));
        for j = [1:i-1, i+1:n]
            l_i = (x - x_node(j)) / (x_node(i) - x_node(j));
            l_x = l_x .* l_i;
        end
        lambda = lambda + abs(l_x);
    end
end

function x_a = define_xa(n)
    x_a = zeros(n+1);
    for i = 0:n
        x_a(i+1) = -1 + (2*i)/n;
    end
end

function x_b = define_xb(n)
    x_b = zeros(n+1);
    for i = 0:n
        x_b(i+1) = cos((2 * i + 1) / (2 * n + 2) * pi);
    end
end

function x_c = define_xc(n)
    x_c = zeros(n+1);
    for i = 0:n
        x_c(i+1) = cos(i * pi / n);
    end
end