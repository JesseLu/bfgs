function [x, err_hist, f_hist] = newton(x0, H, g, f)

tic;
n = length(x0);
x = x0;


% Find step direction.
for k = 1 : 1e2
    dx = -H(x) \ g(x);
    [x, f_hist(k)] = backtrack(f, g, x, dx, 0.1, 0.5);
    err_hist(k) = norm(g(x))^2 / n;
    if (err_hist(k) <= 1e-10)
        break
    end
end


fprintf('Newton Method (%d steps, %1.1e secs) -> funval %e , error %e\n', ...
    length(err_hist), toc, f_hist(end), err_hist(end));
subplot 211; semilogy(f_hist, '.-');
subplot 212; semilogy(err_hist, '.-');
drawnow

