function Uold = MacCormack(x, V, Courant, dh, Uinit, xStartStep, Lstep, xLeft, ...
        xRight, TextSize, t_fin)
    i = 1; % итерационный индекс элементов
    t(i) = 0; % статистика  времени (от итераций)
    dt = Courant * dh / abs(V); % шаг по времени
    plot_interval = round(1 / dt); % интервал вывода графиков
    plot_time = 1; % итерационный индекс графиков
    % Первое приближение:
    Uold = (x >= xStartStep & x <= (xStartStep + Lstep)) * Uinit;
    Unew = Uold;
    Upred = Uold;
    while t(i) <= t_fin
        % Время итерации:
        t(i+1) = t(i) + dt;
        % Точное решение:
        Uexac = Exact(x, t(i+1), V, Uinit, xStartStep, Lstep, xLeft, xRight);
        % Схема Мак-Кормака (двухшаговая):
        % Первый шаг - предиктор (находится оценка величины U на (n+1)-м шаге
        % по времени, используются разности назад):
        Upred(2:end-1) = Uold(2:end-1) - V * dt / dh * (Uold(2:end-1) - Uold(1:end-2));
        % Второй шаг - корректор (определяется окончательное значение U на (n+1)-м шаге
        % по времени, используются разности вперед):
        Unew(2:end-1) = 0.5 * (Uold(2:end-1) + Upred(2:end-1)) - ...
            V * dt / (2 * dh) * (Upred(3:end) - Upred(2:end-1));
        % Периодические граничные условия:
        Upred(1) = Uold(1) - V * dt / dh * (Uold(1) - Uold(end));

Upred(end) = Uold(end) - V * dt / dh * (Uold(end) - Uold(end-1));
        Unew(1) = 0.5 * (Uold(1) + Upred(1)) - V * dt / dh * (Upred(2) - Upred(1));
        Unew(end) = 0.5 * (Uold(end) + Upred(end)) - ...
            V * dt / dh * (Upred(1) - Upred(end));
        % Конец тела цикла:
Uold = Unew;
        % График изменения решения:
        if i == plot_time
            plot_time = plot_time + plot_interval;
            figure(1), clf;
            plot(x, Uold, '-r*', x, Uexac, '-k', 'LineWidth', 3, 'MarkerSize', 5);
            grid on;
            title(['Mac-Cormack, ', 'Time = ', num2str(t(i+1)), ', dt = ', ...
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