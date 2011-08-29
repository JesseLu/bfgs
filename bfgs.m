function [x, err_hist, f_hist] = bfgs(x0, g, f)
tic

n = length(x0);
x = x0;


% Initial approximate inverse Hessian.
H = eye(n);

for k = 1 : 200
    dx = -H * g(x);
    [x_next, f_hist(k)] = backtrack(f, g, x, dx, 0.1, 0.5);

    err_hist(k) = norm(g(x))^2 / n;
    if (err_hist(k) <= 1e-10)
        break
    end
    s = x_next - x;
    y = g(x_next) - g(x);

    p = 1 / (y' * s);
    if p <= 0
        error('%e', p);
    end

    V = eye(n) - p * y * s';
    H = V' * H * V + p * s * s';

    x = x_next;
end

fprintf('BFGS (%d steps, %1.1e secs) -> funval %e , error %e\n', ...
    length(err_hist), toc, f_hist(end), err_hist(end));
subplot 211; semilogy(f_hist, '.-');
subplot 212; semilogy(err_hist, '.-');
drawnow
