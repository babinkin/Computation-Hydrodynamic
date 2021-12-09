clc;
clear;

%Ввод начальных параметров.

TextSize = 15;
a = 1; % коэффициент переноса
N = 301; % количество узлов
x_scale = 1; % размер расчетной области
LBV = 0.1; % левое граничное значение
RBV = 0.8; % правое граничное значение
BreakPoint = 0.5; % точка разрыва

dh = x_scale / (N - 1); % шаг по пространству
dt_allowed = dh^2 / (2 * a); % верхняя граница шага по времени


dt_input= inputdlg(['Allowed time step dt <= ', num2str(dt_allowed), ...
        '. Input time step: '], 'Time step', [1 50], {num2str(dt_allowed)});

    dt = str2double(dt_input{1}); % шаг по времени
    sp = 10; % остановочное число
  
    
%Введение сетки и векторов подстановок.

x = 0:dh:x_scale; % сетка по x
LeftPoints = find(x <= BreakPoint*x_scale); % точки слева от точки разрыва
RightPoints = find(x > BreakPoint*x_scale); % точки справа от точки раз


%Явная схема
%Uexpl = Explicit(x, a, N, LBV, RBV, dh, dt, LeftPoints, RightPoints, sp, TextSize);

%Схема Дюфорта-Франкела 
%Uleap = LeapFrogDF(x, a, N, LBV, RBV, dh, dt, LeftPoints, RightPoints, sp, TextSize);

%Схема Аллена-Чена
%Uac = AllenChen(x, a, N, LBV, RBV, dh, dt, LeftPoints, RightPoints, sp, TextSize);

%plot(x, Uexpl, '.b', x, Uleap, '-r', x, Uac, '-g', 'LineWidth', 1);
%grid on;
%title('U(x)', 'FontSize', TextSize);
%xlabel('x', 'FontSize', TextSize);
%ylabel('U', 'FontSize', TextSize);
%legend({'Explicit', 'Leap-Frog (Du Fort — Frankel)', 'Allen — Chen'}, 'Location', 'best');