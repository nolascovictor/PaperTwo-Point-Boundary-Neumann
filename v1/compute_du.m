function du = compute_du(a, x)
  m = length(a);
  ii = 1:(m-1);
  omega = [1 sqrt(1 + (pi*ii).^2 + (pi*ii).^4)];
  du = 0;
  for i = 2:m
    du = du - (((i-1)*pi) * sqrt(2) * sin((i-1) * pi * x) / omega(i)) * a(i);
  end
end
