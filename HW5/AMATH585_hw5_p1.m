



x = [0:0.01:5];
g = @(x)((5 - x) * exp(x) - 5);
fplot(g,[0,5])
grid on
title('f(x) on [0,5]')
ylabel('(5-x)exp(x)-5')
xlabel('x')

a = 4;
b = 5;
tol_b = 1e-6;
[sol_b, an, bn] = bisection(a, b, tol_b);
steps_b = 1:1:size(an);
var_Name = {'Step', 'an', 'bn'};
t_b = table(steps_b.', an, bn, 'VariableNames',var_Name)
sol_b


x0 = 5;
tol_n = 1e-8;
[sol_n, xn, fxn] = newton(x0, tol_n);
steps_n = 1:1:size(xn);
var_Name = {'Step', 'xk', 'f(xk)'};
t_n = table(steps_n.', xn, fxn, 'VariableNames',var_Name)
sol_n

x0 = 4;
x1 = 5;
tol_s = 1e-16;
[sol_s, xn, fxn] = secant(x0, x1, tol_s);
steps_n = 1:1:size(xn);
var_Name = {'Step', 'xk', 'f(xk)'};
t_s = table(steps_n.', xn, fxn, 'VariableNames',var_Name)
sol_s

function y = f(x)
    y = (5 - x) * exp(x) - 5;
end

function dy = df(x)
    dy = (4 - x) * exp(x);
end

function [sol, an, bn] = bisection(a, b, tol)
    an = [a];
    bn = [b];
    while (bn-an)/2 > tol
        xn = (a + b) / 2;
        if f(xn) < 0
            a = xn;
        else
            b = xn;
        end
        an = [an; a];
        bn = [bn; b];
    end
    sol = (bn(end)+an(end))/2;
end

function [sol, xn, fxn] = secant(x0, x1, tol_s)
    xn = [];
    fxn = [];
    while abs(f(x1)) > tol_s
        xn = [xn; x1];
        fxn = [fxn; f(x1)];
        x_holder = x1 - (x1 - x0)/(f(x1)-f(x0))*f(x1);
        x0 = x1;
        x1 = x_holder;
    end
    sol = x1;
end

function [sol, xn, fxn] = newton(xk, tol_n)
    xn = [];
    fxn = [];
    while abs(f(xk)) > tol_n
        xn = [xn; xk];
        fxn = [fxn; f(xk)];
        xk = xk-f(xk)/df(xk);
    end
    sol = xk;
end

% [ntol, x] = sbisec(4, 5, 1e-6)
% function [ntol,x]=sbisec(a,b,tol)
%     ntol=ceil(log((b-a)/tol)/log(2));
%     for n=1:ntol
%         x=(a+b)/2;
%         fx=f(x);
%         disp([a,b])
%         if fx<0
%             a=x;
%         else
%             b=x;
%         end
%     end
% end