 function Uold = TwoStepLaxWendroff(x, V, Courant, dh, Uinit, xStartStep, Lstep, xLeft, ...
        xRight, TextSize, t_fin)
    i = 1; % итерационный индекс элементов
t(i) = 0; % статистика  времени (от итераций)
dt = Courant * dh / abs(V); % шаг по времени
plot_interval = round(1 / dt); % интервал вывода графиков
plot_time = 1; % итерационный индекс графиков
% Первое приближение:
Uold = (x >= xStartStep & x <= (xStartStep + Lstep)) * Uinit;
Unew = Uold;
Uhalf1 = Uold;
Uhalf2 = Uold;
while t(i) <= t_fin
    % Время итерации:
    t(i+1) = t(i) + dt;
    % Точное решение:
    Uexac = Exact(x, t(i+1), V, Uinit, xStartStep, Lstep, xLeft, xRight);
    % Схема Лакса — Вендроффа (двухшаговая):
    % Первый шаг - предиктор (обеспечивает выполнение условий устойчивости):
    Uhalf1(2:end-1) = 0.5 * (Uold(3:end) + Uold(2:end-1)) - ...
        V * dt / (2 * dh) * (Uold(3:end) - Uold(2:end-1));
    Uhalf2(2:end-1) = 0.5 * (Uold(1:end-2) + Uold(2:end-1)) - ...
        V * dt / (2 * dh) * (Uold(2:end-1) - Uold(1:end-2));
    % Второй шаг - корректор (обеспечивает выполнение требуемой точности):
    Unew(2:end-1) = Uold(2:end-1) - V * dt / dh * (Uhalf1(2:end-1) - Uhalf2(2:end-1));
    % Периодические граничные условия:
    Uhalf1(1) = 0.5 * (Uold(2) + Uold(1)) - ...
        V * dt / (2 * dh) * (Uold(2) - Uold(1));
    Uhalf2(1) = 0.5 * (Uold(end) + Uold(1)) - ...
        V * dt / (2 * dh) * (Uold(1) - Uold(end));
    Uhalf1(end) = 0.5 * (Uold(1) + Uold(end)) - ...
        V * dt / (2 * dh) * (Uold(1) - Uold(end));
    Uhalf2(end) = 0.5 * (Uold(end-1) + Uold(end)) - ...
        V * dt / (2 * dh) * (Uold(end) - Uold(end-1));
    Unew(1) = Unew(1) - V * dt / dh * (Uhalf1(1) - Uhalf2(1));
    Unew(end) = Unew(end) - V * dt / dh * (Uhalf1(end) - Uhalf2(end));
    % Конец тела цикла:
Uold = Unew;
    % График изменения решения:
    if i == plot_time
        plot_time = plot_time + plot_interval;
        figure(1), clf;
        plot(x, Uold, '-r*', x, Uexac, '-k', 'LineWidth', 3, 'MarkerSize', 5);
         grid on;
            title(['Lax — Wendroff (2-step), ', 'Time = ', num2str(t(i+1)), ...
                ', dt = ', num2str(dt)], 'FontSize', TextSize);
            xlabel('X');
            ylabel('U');
            xlim([xLeft, xRight]);
            ylim([0, Uinit*1.5]);
            pause(0.05);
end
i = i + 1;
end

 end