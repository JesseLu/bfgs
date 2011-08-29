% Implement backtracking line search.

function [x, f1] = backtrack(f, g, x, dx, alpha, beta)

f0 = f(x); % Initial objective value.
g0 = g(x); % Initial gradient.
t = 1; % Initial step size.

for k = 1 : 100
    f1 = f(x + t * dx);
    if (f1  <= (f0 - alpha * t * real(g0' * dx)))
        % Step length should satisfy Wolfe conditions.
        x = x + t * dx;
        return
    else
        t = beta * t;
    end
end

error('Backtracking line search failed.');
