function [Uold, q_l, t] = CrankNicolson1(x, lambda, c, a, alpha, T_init, T_gaz, delta, t_fin, N, dh, dt, TextSize)
    
    % Начальные условия:
    Uold = T_init * ones(1, N);
    i = 1;
    t(i) = 0;
    q_l(i) = abs(alpha * (T_gaz - Uold(1)));
    
    % Коэффициенты разностной схемы:
    Bi = alpha * dh / lambda; % кооэффициент Био (безразмерный коэффициент теплоотдачи)
    LBV = Bi * T_gaz / (1 + Bi); % левое граничное значение
    RBV = 0; % правое граничное значение
    A = [ones(1, N-2), -1];
    B = [1, -2*(1+dh^2/(a*dt))*ones(1, N-2), 1];
    C = [-1/(1+Bi), ones(1, N-2)];
    M = full(gallery('tridiag', A, B, C));
   
    while t(i) <= t_fin
       
        t(i+1) = t(i) + dt;
        D = ([LBV, -1*Uold(1:N-2)+2*(1-dh^2/(a*dt))*Uold(2:N-1)-1*Uold(3:N), RBV])';
        Unew =  mldivide(M, D);
        Uold = Unew';
       
        % Вывод графиков:
        plot(x, Unew, 'b', 'LineWidth', 1);
        grid on;
        title(['Crank — Nicolson (1): ', newline, '\delta = ', num2str(delta), ...
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
