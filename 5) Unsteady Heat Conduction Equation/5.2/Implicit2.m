function [Uold, q_l, t] = Implicit2(x, lambda, c, a, alpha, T_init, T_gaz, delta, t_fin, N, dh, dt, TextSize)
    
    % Начальные условия:
    Uold = T_init * ones(1, N);
    i = 1;
    t(i) = 0;
    q_l(i) = abs(alpha * (T_gaz - Uold(1)));
   
    % Коэффициенты разностной схемы:
    Bi = alpha * dh / lambda; % кооэффициент Био (безразмерный коэффициент теплоотдачи)
    LBV = 2 * Bi * T_gaz / (3 + 2 * Bi); % левое граничное значение
    RBV = 0; % правое граничное значение
    A = [ones(1, N-2), -4/3];
    B = [1, -2*(1+dh^2/(2*a*dt))*ones(1, N-2), 1];
    C = [-4/(3+2*Bi), ones(1, N-2)];
    M = full(gallery('tridiag', A, B, C));
    M(1,3) = 1 / (3 + 2 * Bi);
    M(N,N-2) = 1 / 3;
    
    while t(i) <= t_fin
       
        t(i+1) = t(i) + dt;
        D = ([LBV, -(dh^2/(a*dt))*Uold(2:N-1), RBV])';
        Unew = mldivide(M, D);
        Uold = Unew';
        
        % Вывод графиков:
        plot(x, Unew, 'b', 'LineWidth', 1);
        grid on;
        title(['Implicit (2): ', newline, '\delta = ', num2str(delta), ...
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
