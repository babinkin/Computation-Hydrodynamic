function [Uold, q_l, t] = Implicit1(x, lambda, c, a, alpha, T_init, T_gaz, delta, t_fin, N, dh, dt, TextSize)
    
    % Начальные условия:
    Uold = T_init * ones(1, N);
    Unew = Uold;
    i = 1;
    t(i) = 0;
    q_l(i) = abs(alpha * (T_gaz - Uold(1)));
   
    % Коэффициенты разностной схемы:
    A = 1;
    B = 1 + dh^2 / (2 * dt * a);
    C = 1;
    Bi = alpha * dh / lambda; % кооэффициент Био (безразмерный коэффициент теплоотдачи)
    
    while t(i) <= t_fin
        
        % Начальные значения векторов прогоночных коэффициентов:
        K = zeros(1, N-1);
        L = K;
        K(1) = 1 / (1 + Bi);
        L(1) = Bi * T_gaz / (1 + Bi);
        t(i+1) = t(i) + dt;
       
        % Прямая прогонка:
        for j = 2:N-1
           
            K(j) = -1 / (K(j-1) - 2 * B);
            L(j) = (-(dh^2 / (dt * a)) * Uold(j) - L(j-1)) / (K(j-1) - 2 * B);

        end
        
        Unew(N) = L(N-1) / (1 - K(N-1));
        
        % Обратная прогонка
        for j = N-1:-1:1
            
            Unew(j) = K(j) * Unew(j+1) + L(j);

        end
        
        Uold = Unew;
       
        % Вывод графиков:
        plot(x, Unew, 'b', 'LineWidth', 1);
        grid on;
        title(['Implicit (1): ', newline, '\delta = ', num2str(delta), ...
            ', \lambda = ', num2str(lambda), ...
            ', \alpha = ', num2str(alpha), ', c = ', num2str(c)], 'FontSize', TextSize);
        xlabel('X', 'FontSize', TextSize);
        ylabel('U', 'FontSize', TextSize);
        axis([0, delta, T_init, 1.1*T_gaz]);
        pause(0.05);
        
        % Тепловой поток:
        q_l(i+1) = abs(alpha * (T_gaz - Uold(1)));

        i = i + 1;

    end
end