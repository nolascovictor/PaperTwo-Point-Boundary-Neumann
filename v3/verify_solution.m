function [eta, nu, K, t_star] = verify_solution(f, f_u, f_v, f_uu, f_uv, f_vv, b, R)
  % Compute eta
  eta = compute_eta(f, b);

  % Compute nu
  nu = compute_nu(f_u, f_v, b);

  % Compute r and K
  norm_w = compute_norm_w(b);
  r = norm_w;
  K = compute_K(f_uu, f_uv, f_vv, r);

  % Compute t*
  Delta = nu^2 - 2 * eta * K;
  if Delta > 0
    t_star = (nu - sqrt(Delta)) / K;
    if (inf(t_star) > 0) && (sup(t_star) < R)
      disp('Proof was successful!')
      disp(['eta = ', num2str(sup(eta))])
      disp(['nu = ', num2str(inf(nu))])
      disp(['K = ', num2str(sup(K))])
      disp(['t_star = ', num2str(sup(t_star))])
    else
      disp('Proof Failed! Root out of bounds.')
    end
  else
    disp('Proof Failed! No real roots.')
  end
end

function eta = compute_eta(f, b)
  w = @(x) compute_iu(b, x);
  dw = @(x) compute_idu(b, x);
  d2w = @(x) compute_id2u(b, x);
  F_u_2 = @(x) (d2w(x) - f(x, w(x), dw(x))).^2;
  eta = sqrt(abs(int_sympson(F_u_2, 0, 1, 1000)));
  eta = intval(sup(eta));
end

function nu = compute_nu(f_u, f_v, b)
  m = length(b);
  ipi = intval('pi');
  im = intval(m);
  omega_m_p_1 = sqrt(intval('1') + (ipi * im)^2 + (ipi * im)^4);

  nu0 = compute_nu0(f_u, f_v, b);
  N = compute_N(f_u, f_v, b);
  L = min(nu0, (ipi * im)^2 / omega_m_p_1);
  nu = intval(inf(L - (N / (ipi * im))));
end

function nu0 = compute_nu0(df_du, df_ddu, b)
  m = length(b);

  ipi = intval('pi');
  ii = intval(1:(m-1));
  omega = [intval('1') sqrt(1 + (ipi.*ii).^2 + (ipi.*ii).^4)];

  % Compute A = DF_cosm(w)
  w = @(x) compute_iu(b, x);
  dw = @(x) compute_idu(b, x);
  A = intval(zeros(m));

  for j = 1:m
    if j == 1
      uj = @(x) intval(1);
      duj = @(x) intval(0);
      d2uj = @(x) intval(0);
    else
      uj = @(x) sqrt(ii(2)) * cos(ii(j-1) * ipi * x) / omega(j);
      duj = @(x) -(ipi * ii(j-1)) * sqrt(ii(2)) * sin(ii(j-1) * ipi * x) / omega(j);
      d2uj = @(x) -((ipi * ii(j-1))^2) * uj(x);
    end

    for i = 1:m
      if i == 1
        Aij_base = @(x) d2uj(x) - df_du(x, w(x), dw(x)) .* uj(x) - df_ddu(x, w(x), dw(x)) .* duj(x);
      else
        Aij_base = @(x) (d2uj(x) - df_du(x, w(x), dw(x)) .* uj(x) - df_ddu(x, w(x), dw(x)) .* duj(x)) .* sqrt(ii(2)) .* cos(ii(i-1) * ipi * x);
      end
      A(i,j) = int_sympson(Aij_base, 0, 1, 1000);
    end
  end

  % Compute nu
  B = A^(-1);
  % Frobenius norm of B
  nu0 = sqrt(sum(B(:).^2))^(-1);
  nu0 = intval(inf(nu0));
end

function N = compute_N(df_du, df_ddu, b)
  w = @(x) compute_iu(b, x);
  dw = @(x) compute_idu(b, x);
  df_dw = @(x) df_du(x, w(x), dw(x));
  df_ddw = @(x) df_ddu(x, w(x), dw(x));
  % Define a Taylor model in the interval [0, 1]
  % dw{k} is the k-th derivative of the function
  dw = df_dw(taylorinit(infsup(0,1), 1));
  ddw = df_ddw(taylorinit(infsup(0,1), 1));
  if isfloat(dw)
    N_dw = intval(dw);
  else
    N_dw = abs(dw{0}) + sqrt(dw{0}^2 + dw{1}^2);
  end
  if isfloat(ddw)
    N_ddw = intval(ddw);
  else
    N_ddw = abs(ddw{0}) + sqrt(ddw{0}^2 + ddw{1}^2);
  end
  N = intval(sup(N_dw + N_ddw));
end

function norm_w = compute_norm_w(b)
  w = @(x) compute_iu(b, x);
  dw = @(x) compute_idu(b, x);
  d2w = @(x) compute_id2u(b, x);
  w_func = @(x) w(x).^2 + dw(x).^2 + d2w(x).^2;
  norm_w = sqrt(abs(sqrt(int_sympson(w_func, 0, 1, 1000))));
  norm_w = intval(sup(norm_w));
end

function K = compute_K(f_uu, f_uv, f_vv, r)
  ii = intval(1:2);
  c_1 = sqrt((exp(ii(2)) + 1) / (exp(ii(2)) - 1));
  K_uu = f_uu(taylorinit(infsup(inf(-c_1 * r), sup(c_1 * r)), 1));
  K_uv = f_uv(taylorinit(infsup(inf(-c_1 * r), sup(c_1 * r)), 1));
  K_vv = f_vv(taylorinit(infsup(inf(-c_1 * r), sup(c_1 * r)), 1));
  if isfloat(K_uu), K_uu = intval(K_uu); else, K_uu = K_uu{0}; end
  if isfloat(K_uv), K_uv = intval(K_uv); else, K_uv = K_uv{0}; end
  if isfloat(K_vv), K_vv = intval(K_vv); else, K_vv = K_vv{0}; end
  K = intval(sup( K_uu^2 + 2 * K_uv^2 + K_vv^2 ));
end

function X = int_sympson(f, a, b, n)
  D = intval(b) - intval(a); H = D / n;
  x = a + intval(0:n) / n*D;
  w = 2 * ones(1,n+1); w(2:2:n) = 4; w(1) = 1; w(n+1) = 1;
  V = sum(H/3 * w .* f(intval(x)));
  Y = f(taylorinit(infsup(a,b), 4)); % Approximation inclusion
  E = H^4*D / intval('7.5') * Y{4}; % Error term
  X = V - E; % Integral inclusion
end

function iu = compute_iu(a, x)
  m = length(a);
  ipi = intval('pi');
  ii = intval(1:(m-1));
  omega = [intval('1') sqrt(1 + (ipi*ii).^2 + (ipi*ii).^4)];
  iu = a(1);
  for i = 2:m
    iu = iu + (sqrt(ii(2)) * cos((i-1) * ipi * x) / omega(i)) * a(i);
  end
end

function idu = compute_idu(a, x)
  m = length(a);
  ipi = intval('pi');
  ii = intval(1:(m-1));
  omega = [intval('1') sqrt(1 + (ipi*ii).^2 + (ipi*ii).^4)];
  idu = 0;
  for i = 2:m
    idu = idu - ((ii(i-1) * ipi) * sqrt(ii(2)) * sin(ii(i-1) * ipi * x) / omega(i)) * a(i);
  end
end

function id2u = compute_id2u(a,x)
  m = length(a);
  ipi = intval('pi');
  ii = intval(1:(m-1));
  omega = [intval('1') sqrt(1 + (ipi*ii).^2 + (ipi*ii).^4)];
  id2u = 0;
  for i = 2:m
    id2u = id2u - ((ii(i-1) * ipi)^2 * sqrt(ii(2)) * cos(ii(i-1) * ipi * x) / omega(i)) * a(i);
  end
end
