
% Construct h
h_values = [];
estimation = [];
errors = [];
for i = 1:16
    power = 10 ^(-1 * i);
    h_values = [h_values; power];
end
size(h_values)
for h_i = h_values.'
    h_i
    [est, error] = calculation(h_i, pi/6);
    estimation = [estimation; est];
    errors = [errors; error];
end


size(estimation)
size(errors)
var_Name = {'h', 'Computed f"(x)', 'Error'};
t = table(h_values, estimation, errors, 'VariableNames',var_Name)

function [est, error] = calculation(h, x)
    est= (sin(x+h) - 2*sin(x) + sin(x-h))/(h^2)
    error = abs(-0.5 - est)
end
