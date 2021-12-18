function [F1,F2,F3,iteration] = SolveFluxExacRieman(U1L,U2L,U3L,U1R,U2R,U3R,Gamma)

% Переводим в физические переменные:
left_density = U1L;
right_density = U1R;
left_velocity = U2L./left_density;
right_velocity = U2R./right_density;
left_pressure = (Gamma - 1) * (U3L - 0.5 * (left_velocity.^2).*left_density);
right_pressure = (Gamma - 1) * (U3R - 0.5 * (right_velocity.^2).*right_density);

MaxIteration = 20; % максимальное число итераций
TOL = 1e-8;
lambda = 0; % линия на грани КО

len = length(U1L);
rho = zeros(1, len);
u = zeros(1, len);
p = zeros(1, len);
for i = 1:len
    [rho(i), u(i), p(i), ~, iteration] = ExacRiemanSolver(left_density(i), ...
        left_velocity(i), left_pressure(i), right_density(i), right_velocity(i), ...
        right_pressure(i), Gamma, lambda, MaxIteration, TOL);

end

% Рассчитываем потоки:
F1 = rho.*u;
F2 = (rho.*u).*u + p;
F3 = (p./(Gamma - 1) + rho.*(u.^2) / 2 + p).*u;

end
% Входные параметры:
% U1L, U2L, U3L - вектор консервативных переменных слева,
% U1R, U2R, U3R - вектор консервативных переменных справа.
%
% Выходные параметры:
% F1, F2, F3 - вектор потока на грани КО,
% iteration - выполненое количество итераций.