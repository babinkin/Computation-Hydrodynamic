function Uold = LaxWendroff(x, V, Courant, dh, Uinit, xStartStep, Lstep, xLeft, ...
        xRight, TextSize, t_fin)
    i = 1; % итерационный индекс элементов
    t(i) = 0; % статистика  времени (от итераций)
    dt = Courant * dh / abs(V); % шаг по времени
    plot_interval = round(1 / dt); % интервал вывода графиков
    plot_time = 1; % итерационный индекс графиков

% Первое приближение:
    Uold = (x >= xStartStep & x <= (xStartStep + Lstep)) * Uinit;
    Unew = Uold;
    while t(i) <= t_fin
        % Время итерации:
        t(i+1) = t(i) + dt;
        % Точное решение:
        Uexac = Exact(x, t(i+1), V, Uinit, xStartStep, Lstep, xLeft, xRight);
        % Схема Лакса — Вендроффа:
        Unew(2:end-1) = Uold(2:end-1) - ...
            V * dt / (2 * dh) * (Uold(3:end) - Uold(1:end-2)) + ...
            V^2 * dt^2 / (2 * dh^2) * (Uold(3:end) - 2 * Uold(2:end-1) + Uold(1:end-2));
        % Периодические граничные условия:
        Unew(1) = Unew(1) - V * dt / (2 * dh) * (Uold(2) - Uold(end)) + ...
            V^2 * dt^2 / (2 * dh^2) * (Uold(2) - 2 * Uold(1) + Uold(end));
        Unew(end) = Unew(end) - V * dt / (2 * dh) * (Uold(1) - Uold(end-1)) + ...
            V^2 * dt^2 / (2 * dh^2) * (Uold(1) - 2 * Uold(end) + Uold(end-1));
        % Конец тела цикла:
Uold = Unew;
        % График изменения решения:
        if i == plot_time
            plot_time = plot_time + plot_interval;
            figure(1), clf;
            plot(x, Uold, '-r*', x, Uexac, '-k', 'LineWidth', 3, 'MarkerSize', 5);
            grid on;
            title(['Lax — Wendroff, ', 'Time = ', num2str(t(i+1)), ', dt = ', ...
                num2str(dt)], 'FontSize', TextSize);
            xlabel('X');
            ylabel('U');
            xlim([xLeft, xRight]);
            ylim([0, Uinit*1.5]);
            pause(0.05);
end
i = i + 1;
    end
end