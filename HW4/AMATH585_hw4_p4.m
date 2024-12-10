% AMATH 585 HW4 Problem 4
% Tianbo Zhang 1938501
% Implementation for Romberg Extrapolation in probolem 4

% Set up the conditions
a = 0;
b = 1;
n = 100;
tol = 1e-12;
[integral, num_eval] = rombergIntegration(a, b, n, tol);
fprintf('Integral: %.15f\n', integral);
fprintf('Number of function evaluations: %d\n', num_eval);

function [integral, num_eval] = rombergIntegration(a,b,n, tol)
% Set up the base condition
    num_eval = 2;
    T=zeros(n,n);
    h=b-a;
    T(1,1)=h*(f(a)+f(b))/2;
    m = 1;
    for i=2:n
        % Trapezoidal rule
        h=h/2;
        sum = 0;
        m = 2*m;
        mm = m-1;
        for k = (1:2:mm)
            sum = sum + f(a+k*h);
        end
        T(i,1)=T(i-1,1)/2+h*sum;

        % Romberg extrapolation
        for k=2:i
            T(i,k)=T(i,k-1)+(T(i,k-1)-T(i-1,k-1))/(4^(k-1)-1);
        end

        % Check for tolerance
        num_eval = num_eval + 2^(i-2);
        integral = T(i,i);
        if abs(T(i,i) - T(i-1, i-1)) < tol
            return;
        end
    end
end

% Function to integrate
function y = f(x)
    y = cos(x.^2);
end