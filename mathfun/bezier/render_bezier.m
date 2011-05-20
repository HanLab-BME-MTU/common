function Y = render_bezier(P, t)

n = size(P, 1) - 1;
B = bsxfun(@power, t, 0:n) .* bsxfun(@power, 1 - t, n:-1:0);
B = B * diag([1 cumprod(n:-1:1) ./ cumprod(1:n)]);
Y = B * P;
