function u = compute_u(a, x)
  m = length(a);
  ii = 1:(m-1);
  omega = [1 sqrt(1 + (pi*ii).^2 + (pi*ii).^4)];
  u = a(1);
  for i = 2:m
    u = u + sqrt(2).*cos((i-1).*pi.*x) ./ omega(i).*a(i);
  end
end
