function [b, fval] = compute_solution(f, b0)
  err_tol = 1e-14;
  [b, fval] = fsolve(@(b) F_cosm(f, b), b0);
  b(abs(b) < err_tol) = 0;
end

function F_a = F_cosm(f, a)
  m = length(a);
  F_u = @(x) inf(compute_d2u(a,x) - f(x, compute_u(a,x), compute_du(a,x)));
  F_a = zeros(1, m);
  F_a(1) = integral(@(x) F_u(x), 0, 1);
  for i = 2:m
    F_a(i) = integral(@(x) F_u(x) * sqrt(2) .* cos((i-1) * pi * x), 0, 1);
  end
end

function d2u = compute_d2u(a, x)
  m = length(a);
  ii = 1:(m-1);
  omega = [1 sqrt(1 + (pi*ii).^2 + (pi*ii).^4)];
  d2u = 0;
  for i = 2:m
    d2u = d2u - (((i-1).*pi).^2 .* sqrt(2) * cos((i-1).*pi.*x) ./ omega(i)) .* a(i);
  end
end
