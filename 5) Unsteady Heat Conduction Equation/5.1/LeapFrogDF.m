function Uold2 = LeapFrogDF(x, a, N, LBV, RBV, dh, dt, LeftPoints, RightPoints, sp, TextSize)
   
    t = 0; % переменная времени
    plot_interval = 1000; % интервал вывода графиков
    plot_time = 1; % итерационный индекс графиков
    step = 0; % итерационный индекс цикла
    stop = sp * plot_interval; % остановочное значение
    tmp1 = 2 * dt * a / dh^2; % вспомогательная переменная счета
    tmp2 = 1 + tmp1; % вспомогательная переменная счета
    
    % Начальные значения:
    Uold1 = zeros(1, N);
    Uold2 = Uold1;
    Unew = Uold1;
 
    Uold1(LeftPoints) = LBV;
    Uold1(RightPoints) = RBV;
    Uold2(LeftPoints) = LBV;
    Uold2(RightPoints) = RBV;
    
    while 1
       
        step = step + 1;
        t = t + dt;
       
        % График изменения решения:
        if step == plot_time
            
            plot_time = plot_time + plot_interval;
            plot(x, Uold2, 'r', 'LineWidth', 2);
            grid on;
            title(['Leap-Frog (Du Fort — Frankel), t = ', num2str(t)], ...
                'FontSize', TextSize);
            xlabel('x', 'FontSize', TextSize);
            ylabel('U', 'FontSize', TextSize);
            pause(0.5);

        end
        
        % Трехслойная схема Дюфорта — Франкела (схема "чехарда"):
        Unew(2:N-1) = tmp1 / tmp2 * (Uold2(1:N-2) - Uold1(2:N-1) + Uold2(3:N)) + Uold1(2:N-1) / tmp2;
        
        % Граничные условия:
        Unew(1) = LBV;
        Unew(N) = RBV;
       
        if step >= stop
            
            break;

        end
        
        Uold1 = Uold2;
        Uold2 = Unew;

    end
end
