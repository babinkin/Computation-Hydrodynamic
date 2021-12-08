clc;
clear;
TextSize = 15;

%нагрев тела из метала

%lambda = 30; % теплопроводность материала стенки, Вт/(м*К)
%c = 300; % теплоемкость материала стенки, Дж/(кг*К)
%rho = 7800; % плотность материала стенки, кг/м^3

% Малотеплопроводный материал

lambda = 1; % теплопроводность материала стенки, Вт/(м*К)
c = 900; % теплоемкость материала стенки, Дж/(кг*К)
rho = 1600; % плотность материала стенки, кг/м^3

a = lambda / (rho * c); % коэффициент температуропроводности
alpha = 4000; % коэффициент теплоотдачи, Вт/(м^2*К)
T_init = 300; % начальная температура стенки, К
T_gaz = 2600; % температура потока, К
delta = 0.005; % толщина стенки, м
t_fin = 2.5; % продолжительность нагрева, с
time_steps = 50; % количество шагов по времени
N = 200; % количество узлов сетки
dh = delta / (N - 1); % шаг по пространству
dt = t_fin / time_steps; % шаг по времени
theta = 0.5;

%Назначение узлов сетки.
x = 0:dh:delta; % сетка по x

%Неявная схема с 1 аппроксимацией граничных условий
[Uimpl1, q_impl1, t_i1] = Implicit1(x, lambda, c, a, alpha, T_init, T_gaz, delta, t_fin, N, dh, dt, TextSize);

%Неявная схема со 2 аппроксимацией граничных условий
[Uimpl2, q_impl2, t_i2] = Implicit2(x, lambda, c, a, alpha, T_init, T_gaz, delta, t_fin, N, dh, dt, TextSize);

%Схема Кранка-Николсона (1 аппроксимация):
[Ucn1, q_cn1, t_cn1] = CrankNicolson1(x, lambda, c, a, alpha, T_init, T_gaz, delta, t_fin, N, dh, dt, TextSize);

%Схема Кранка-Николсона (2 аппроксимация):
[Ucn2, q_cn2, t_cn2] = CrankNicolson2(x, lambda, c, a, alpha, T_init, T_gaz, delta, t_fin, N, dh, dt, TextSize);

%Схема Рихтмайера-Мортона (1 аппроксимация):
[Urm1, q_rm1, t_rm1] = RichtmyerMorton1(x, lambda, c, a, alpha, T_init, T_gaz, delta, t_fin, N, dh, dt, theta, TextSize);

%Схема Рихтмайера-Мортона (2 аппроксимация):
[Urm2, q_rm2, t_rm2] = RichtmyerMorton2(x, lambda, c, a, alpha, T_init, T_gaz, delta, t_fin, N, dh, dt, theta, TextSize);

plot(t_i1, q_impl1, t_i2, q_impl2, t_cn1, q_cn1, t_cn2, q_cn2, ...
    t_rm1, q_rm1, t_rm2, q_rm2,  'LineWidth', 1);
grid on;
title('Heat flow', 'FontSize', TextSize);
xlabel('t, с', 'FontSize', TextSize);
ylabel('q_l, Вт/м^2', 'FontSize', TextSize);
axis([0, t_fin, 0, 10e6]);
legend({'Implicit (1)', 'Implicit (2)', 'Crank — Nicolson (1)', ...
    'Crank — Nicolson (2)', 'Richtmyer — Morton (1)', ...
    'Richtmyer — Morton (2)'}, 'Location', 'best');
