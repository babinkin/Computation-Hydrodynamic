clc
close all
clear

TextSize = 20;

% Вводим начальные данные
V = 1.005; % скорость переноса
Uinit =  1; % высота начальной ступеньки
xStartStep = 0.0; % координата от которой начинается ступенька
Lstep = 5.0; % ширина начальной ступеньки
xLeft = 0; % левая граница расчетной области
xRight = 100; % правая граница расчетной области
% ----
Courant = 0.9; % число Куранта
N = 101; % число узлов сетки
dh = (xRight-xLeft)/(N-1); % шаг сетки
t_finish = 50.0;
% --------
% задаем сетку
x = xLeft:dh:xRight;
% ------
% Начальные условия - ступенька
%Uold = (x>=xStartStep & x<=(xStartStep + Lstep))*Uinit;
% Uold = U_Exac_Step(x,0,V,Uinit,xStartStep,Lstep,xLeft,xRight);
% ---
% выводим первое приближение
%plot(x,Uold,'-k*','LineWidth',3,'MarkerSize',5)
%grid on
%xlim([xLeft xRight])
%ylim([0 Uinit*1.5])
%xlabel('X')
%ylabel('U')
%set(gca,'FontSize',20)
%pause(0.5)
% return 
% ---------------------------------


%U_left = LeftAngle(x, V, Courant, dh, Uinit, xStartStep, Lstep, xLeft, xRight, ...
    %TextSize, t_finish);


%U_towards = Towards(x, V, Courant, dh, Uinit, xStartStep, Lstep, xLeft, xRight, ...
   % TextSize, t_finish);

%U_central = Central(x, V, Courant, dh, Uinit, xStartStep, Lstep, xLeft, xRight, ...
    %TextSize, t_finish);

%U_Lax = Lax(x, V, Courant, dh, Uinit, xStartStep, Lstep, xLeft, xRight, ...
   % TextSize, t_finish);

U_MacCormack= MacCormack(x, V, Courant, dh, Uinit, xStartStep, Lstep, xLeft, xRight, ...
    TextSize, t_finish);








































% Основной цикл
NStep = 90; % число шагов
plot_interval = 1; % интервал вывода графиков в шагах
t_fin = 1000.0; % Если досдигнуто время выхода то остановка

% ----
t(1) = 0; % статистика  времени (от итераций)
plot_time=plot_interval;
Unew = Uold;

for i=1:NStep
    i
    dt = Courant*dh/abs(V); % шаг по времени
    t(i+1) = t(i) + dt; 
    % точное решение
    Uexac = U_Exac_Step(x,t(i+1),V,Uinit,xStartStep,Lstep,xLeft,xRight);
    % Реализует шаг по времени "Левый уголок"
    Unew(2:end-1) = Uold(2:end-1) - V*dt/dh*(Uold(2:end-1) - Uold(1:end-2));
    % Схема Лакса
    % ----
    % Граничные условия задаются как переодичные ГУ
    Unew(1) = Unew(1) - V*dt/dh*(Uold(1) - Uold(end));
    Unew(end) = Unew(end) - V*dt/dh*(Uold(end) - Uold(end-1));
    % Схема Лакса ГУ
    % -----
    % Конец тела цикла
    Uold = Unew;
    % ----
    if i==plot_time
        figure(1), clf
        plot_time=plot_time+plot_interval;
        plot(x,Uold,'-r*',x,Uexac,'-k','LineWidth',3,'MarkerSize',5)
        grid on
        xlim([xLeft xRight])
        ylim([0 Uinit*1.5])
        xlabel('X')
        ylabel('U')
        title(['Time = ' num2str(t(i+1)) ', dt=' num2str(dt) ' '],'FontSize',TextSize)
        set(gca,'FontSize',20)
        pause(0.05)
    end
    % ---
    if t(i+1)>t_fin
            break
    end
    
end






































