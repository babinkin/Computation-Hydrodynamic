function [F1,F2,F3,iteration] = SolveFluxExacRieman(U1L,U2L,U3L,U1R,U2R,U3R,Gamma)

% ��������� � ���������� ����������:
left_density = U1L;
right_density = U1R;
left_velocity = U2L./left_density;
right_velocity = U2R./right_density;
left_pressure = (Gamma - 1) * (U3L - 0.5 * (left_velocity.^2).*left_density);
right_pressure = (Gamma - 1) * (U3R - 0.5 * (right_velocity.^2).*right_density);

MaxIteration = 20; % ������������ ����� ��������
TOL = 1e-8;
lambda = 0; % ����� �� ����� ��

len = length(U1L);
rho = zeros(1, len);
u = zeros(1, len);
p = zeros(1, len);
for i = 1:len
    [rho(i), u(i), p(i), ~, iteration] = ExacRiemanSolver(left_density(i), ...
        left_velocity(i), left_pressure(i), right_density(i), right_velocity(i), ...
        right_pressure(i), Gamma, lambda, MaxIteration, TOL);

end

% ������������ ������:
F1 = rho.*u;
F2 = (rho.*u).*u + p;
F3 = (p./(Gamma - 1) + rho.*(u.^2) / 2 + p).*u;

end
% ������� ���������:
% U1L, U2L, U3L - ������ �������������� ���������� �����,
% U1R, U2R, U3R - ������ �������������� ���������� ������.
%
% �������� ���������:
% F1, F2, F3 - ������ ������ �� ����� ��,
% iteration - ���������� ���������� ��������.