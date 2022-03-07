% The function verify_solution works for any f such that
% the second derivatives depend on only one variable

set(0, 'DefaultAxesFontSize', 18)
set(0, 'DefaultAxesFontWeight', 'bold')

% clear; clc;
tic

m = 5; % Projection dimension

% Define initial guess
% b0 = 10 * rand(1, m);
b0 = [0, 0.5, 0, 0, 0];

% Define the equation
ipi = intval('pi');
i2 = intval(2);
i6 = intval(6);
% f(x, u, v), where u = u(x) and v(x) = u'(x)
f = @(x, u, v) -u + (u.^3) / i6 - cos(ipi*x);
f_u = @(x, u, v) -1 + (u.^2) / i2;
f_v = @(x, u, v) 0;
f_uu = @(u) u;
f_uv = @(u) 0;
f_vv = @(u) 0;

% Compute numerical solution
b = compute_solution(f, b0);

% Verify the numerical solution
R = 1;
[eta, nu, K, t_star, t_double_star] = verify_solution(f, f_u, f_v, f_uu, f_uv, f_vv, b, R);

% Plot the solution
w = @(x) compute_u(b, x);
x = 0:0.01:1;
plot(x, w(x), 'r', 'LineWidth', 2);

% axis([0 1 -0.15 0.16])

toc
