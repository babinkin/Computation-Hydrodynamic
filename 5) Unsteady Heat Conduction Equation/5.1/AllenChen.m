function Uold = AllenChen(x, a, N, LBV, RBV, dh, dt, LeftPoints, RightPoints, sp, TextSize)
    
    t = 0; % переменная времени
    plot_interval = 1000; % интервал вывода графиков
    plot_time = 1; % итерационный индекс графиков
    step = 0; % итерационный индекс цикла

    stop = sp * plot_interval; % остановочное значение
    tmp1 = dt * a / dh^2; % вспомогательная переменная счета
    tmp2 = 1 + 2 * tmp1; % вспомогательная переменная счета

    % Начальные значения:
    Uold = zeros(1, N);
    Upred = Uold;
    Unew = Uold;
    Uold(LeftPoints) = LBV;
    Uold(RightPoints) = RBV;
    Upred(1) = LBV;
    Upred(N) = RBV;

    while 1
    step = step + 1;
    t = t + dt;
   
    % График изменения решения:
    if step == plot_time
        
        plot_time = plot_time + plot_interval;
        plot(x, Uold, '-g', 'LineWidth', 2);
        grid on;
        title(['Allen — Chen, t = ', num2str(t)], 'FontSize', TextSize);
        xlabel('x', 'FontSize', TextSize);
        ylabel('U', 'FontSize', TextSize);
        pause(0.5);

    end
    
    % Схема Аллена — Чена (двухшаговая):
    % Первый шаг - предиктор:
    Upred(2:N-1) = Uold(2:N-1) / tmp2 + tmp1 / tmp2 * (Uold(1:N-2) + Uold(3:N));
    
    % Второй шаг - корректор:
    Unew(2:N-1) = Uold(2:N-1) / tmp2 + tmp1 / tmp2 * (Upred(1:N-2) + Upred(3:N));
    
    % Граничные условия:
    Unew(1) = LBV;
    Unew(N) = RBV;
   
    if step >= stop
       
        break;

    end
    
    Uold = Unew;
    
    end
end