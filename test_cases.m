function test_cases(n)

opt{1}.name = 'Newton';
opt{1}.fun = @(x0, H, g, f) newton(x0, H, g, f);

opt{2}.name = 'BFGS';
opt{2}.fun = @(x0, H, g, f) bfgs(x0, g, f);

for k = 1 : length(n)
    fprintf('\nn = %d\n===\n', n(k));
    for l = 2 : length(opt)
        quad(opt{l}, n(k))
        my_log(opt{l}, n(k))
    end
end

function [z] = cust_log(x)
z = log(x);
ind = find(x <= 0);
z(ind) = -Inf;


function my_log(opt, n)


% Create problem.
m = n + randi(n);
A = randn(m, n);
b = randn(m, 1);

f = @(x) -sum(cust_log(x)) + 0.5 * norm(A*x - b)^2;
g = @(x) -x.^-1 + A'*(A*x - b);
H = @(x) diag(x.^-2) + A'*A;

% Pick random starting point.
x0 = abs(randn(n, 1));

% Optimize.
fprintf('Running MY_LOG: ')
[x, err_hist] = opt.fun(x0, H, g, f);


function quad(opt, n)

% Create problem.
A = randn(n) + 1 * randn(n);
b = randn(n, 1);

H = @(x) A' * A; % Hessian.
g = @(x) A' * (A * x - b); % Gradient.
f = @(x) 0.5 * norm(A * x - b)^2; % Objective function.

% Pick random starting point.
x0 = randn(n, 1);

% Optimize.
fprintf('Running QUAD: ')
[x, f_hist] = opt.fun(x0, H, g, f);

