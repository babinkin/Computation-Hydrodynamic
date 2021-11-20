clc
close all
clear
%==========================================================================
% Настройка шрифтов в зависимости от ОС
if isunix
    fontname = 'Free Helvetian';% for LINUX
    % Для Юникс системвыставляется не тот шрифт в функции msgbox, поэтому
    % либо правим файл toolbox/matlab/uitools/msgbox.m, либо для каждого
    % объекта после создания меняем шрифт
elseif ispc
    fontname = 'Arian Cyr';% for Windows
end
set(0,'DefaultAxesFontName',fontname);
set(0,'DefaultTextFontName',fontname);
set(0,'DefaultUIControlFontname',fontname);
set(0,'fixedwidthfontname',fontname);
TextSize = 20;
%==========================================================================

% Входные параметры задачи
h0 = 80000; % высота точки начала интегрирования (точки входа)
theta0 = -4*pi/180; % угол входа
v0 = 7800; % скорость входа
g0 = 9.81; % ускорение свободного падения на поверхности земли
rho0 = 1.23; % плотность воздуха на поверхности земли
K = 0.3; % а/д качесто
rz = 6.37e6; % радиус земли
sigma = 0.01; % баллистический параметр Cx*S/m


% выбираем шаг интегрирования
dt = 5;
% начальные условия
t0 = 0;
x0=0;
% - вектор входных параметров
Y0 = [v0 theta0 x0 h0];


% Выполняем интегрирование системы уравнений, пока не достигнед земли h=0
i = 1; % индекс элемента
Y1(i,:) = Y0;
T1(i) = t0;
while Y1(i,4) > 0
    T1(i+1) = T1(i) + dt;
    Y1(i+1,:) = Y1(i,:) + dt * rightPart(T1(i), Y1(i,:), rho0, g0, rz, sigma, K)';
    i = i + 1;
end




% метод Эйлера второго порядка
 i = 1; % индекс элемента
Y2(i,:) = [v0 theta0 x0 h0]; % вектор входных параметров
T2(i) = t0;
while Y2(i,4) > 0
    T2(i+1) = T2(i) + dt;
    Y2(i+1,:) = Y2(i,:) + dt * (rightPart(T2(i), Y2(i,:), rho0, g0, rz, sigma, K)' + rightPart(T2(i), Y2(i,:) + dt * rightPart(T2(i), Y2(i,:), rho0, g0, rz, sigma, K)', rho0, g0, rz, sigma, K)') / 2;
    i = i + 1;
end


% Решение Библиотечной функцией
[T, Y] = ode45(@(t, y) rightPart(t, y, rho0, g0, rz, sigma, K), linspace(t0, max(T2), length(T2)), Y0);


figure(1);
plot(T1, Y1(:,4)/1000, '-b', T2, Y2(:,4)/1000, 'ro', T, Y(:,4)/1000, '-g');
grid on;
title('h(t)', 'FontSize', TextSize);
xlabel('Время t, с', 'FontSize', TextSize);
ylabel('Высота h, км','FontSize', TextSize);
xlim([0, max(T2)]);
ylim([0, 1.5*h0/1000]);
legend({'Метод Эйлера 1-го порядка', 'Метод Эйлера 2-го порядка', 'ode45'}, ...
    'Location', 'best');


figure(2);
plot(Y1(:,3)/1000, Y1(:,4)/1000, '-b', Y2(:,3)/1000, Y2(:,4)/1000, 'ro', ...
    Y(:,3)/1000, Y(:,4)/1000, '-g');
grid on;
title('h(x)', 'FontSize', TextSize);
xlabel('Путь x, км', 'FontSize', TextSize);
ylabel('Высота h, км', 'FontSize', TextSize);
ylim([0, 1.5*h0/1000]);
legend({'Метод Эйлера 1-го порядка', 'Метод Эйлера 2-го порядка', 'ode45'}, ...
    'Location', 'best');

figure(3);
plot(T1, Y1(:,3)/1000, '-b', T2, Y2(:,3)/1000, 'ro', T, Y(:,3)/1000, '-g');
grid on;
title('x(t)', 'FontSize', TextSize);
xlabel('Время t, с', 'FontSize', TextSize);
ylabel('Путь x, км', 'FontSize', TextSize);
xlim([0, max(T2)]);
legend({'Метод Эйлера 1-го порядка', 'Метод Эйлера 2-го порядка', 'ode45'}, ...
    'Location', 'best');

figure(4);
plot(T1, Y1(:,2)*180/pi, '-b', T2, Y2(:,2)*180/pi, 'ro', T, Y(:,2)*180/pi, '-g');
grid on;
title('\theta(t)', 'FontSize', TextSize);
xlabel('Время t, с', 'FontSize', TextSize);
ylabel('Угол наклона траектории \theta, град.', 'FontSize', TextSize);
xlim([0, max(T2)]);
legend({'Метод Эйлера 1-го порядка', 'Метод Эйлера 2-го порядка', 'ode45'}, ...
    'Location', 'best');

figure(5);
plot(T1, Y1(:,1), '-b', T2, Y2(:,1), 'ro', T, Y(:,1), '-g');
grid on;
title('V(t)', 'FontSize', TextSize);
xlabel('Время t, с', 'FontSize', TextSize);
ylabel('Скорость V, м/с', 'FontSize', TextSize);
xlim([0, max(T2)]);
legend({'Метод Эйлера 1-го порядка', 'Метод Эйлера 2-го порядка', 'ode45'}, ...
    'Location', 'best');



j_left = -10; % левая граница варьирования
j_right = 2; % правая граница варьирования
q = zeros(500, abs(j_right-j_left)); % скоростной напор (кг / (м * с^2)
n_max = zeros(abs(j_right-j_left), 1); % максимальная перегрузка
c = 1; % индекс элемента
for j = j_left:j_right
    theta0 = j * pi / 180; % угол входа, рад
    Y0 = [v0 theta0 x0 h0]; % вектор входных параметров
    [TT, YY] = ode45(@(t, y) rightPart(t, y, rho0, g0, rz, sigma, K), ...
        linspace(t0, 3000, 500), Y0);
    q(:,c) = (0.5 * rho0).*exp(-YY(:,4)./7800).*(YY(:,1).^2);
    n = sqrt(((sigma.*q(:,c))).^2 + (((sigma * K).*q(:,c))).^2) / g0; % перегрузка
    n_max(c) = max(n);
if j == 0
        figure(5);
        plot(TT, n, '-b' );

   grid on;
        title('n(t)', 'FontSize', TextSize);
        xlabel('Время t, с', 'FontSize', TextSize);
        ylabel('Перегрузка n', 'FontSize', TextSize);
        legend({'\theta = 0', 'Max'}, 'Location', 'best');
        xlim([0, max(T2)]);
end
c = c + 1;
end

[n_min, n_min_ind] = min(n_max);

fprintf('При θ = %d°  минимальная перегрузка:  %.2f g.', ...
    j_left+n_min_ind-1, n_min);



















