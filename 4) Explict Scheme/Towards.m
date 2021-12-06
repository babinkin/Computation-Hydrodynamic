function Uold = Towards(x, V, Courant, dh, Uinit, xStartStep, Lstep, xLeft, ...
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
        Uexac = Exact(x, t((V+abs(V))/(2 * V)*(i+1)+(V-abs(V))/(2 * V)*i), V, ...
            Uinit, xStartStep, Lstep, xLeft, xRight);

 % Схема "навстречу потоку":
        Unew(2:end-1) = Uold(2:end-1) - ...
            V * dt / (2 * dh) * (Uold(3:end) - Uold(1:end-2)) + ...
            abs(V) * dt / (2 * dh) * (Uold(1:end-2) - 2 * Uold(2:end-1) + Uold(3:end));
        % Периодические граничные условия:
        Unew(1) = Unew(1) - V * dt / (2 * dh) * (Uold(2) - Uold(end)) + ...
            abs(V) * dt / (2 * dh) * (Uold(end) - 2 * Uold(1) + Uold(2));
        Unew(end) = Unew(end) - V * dt / (2 * dh) * (Uold(1) - Uold(end-1)) + ...
            abs(V) * dt / (2 * dh) * (Uold(end-1) - 2 * Uold(end) + Uold(1));
        % Конец тела цикла:
Uold = Unew;
        % График изменения решения:
        if i == plot_time
            plot_time = plot_time + plot_interval;
            figure(1), clf;
            plot(x, Uold, '-r*', x, Uexac, '-k', 'LineWidth', 3, 'MarkerSize', 5);
            grid on;
            title(['Upstream, ', 'Time = ', ...
                num2str(t((V+abs(V))/(2 * V)*(i+1)+(V-abs(V))/(2 * V)*i)), ...
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