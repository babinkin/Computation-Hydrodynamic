clc
close
clear
if isunix
    fontname = 'Free Helvetian';
  elseif ispc
    fontname = 'Arian Cyr';
end
set(0,'DefaultAxesFontName',fontname);
set(0,'DefaultTextFontName',fontname);
set(0,'DefaultUIControlFontname',fontname);
set(0,'fixedwidthfontname',fontname);
TextSize = 20;



% Количество узлов 
N = 100;
% Диапазон изменения x
x_left = 0;
x_right = 5;
% вычисляем шаг по x (равномерный)
h = (x_right - x_left)/(N - 1);

% получаем значение функции в узловых точках
x = x_left:h:x_right;
y = f_ot_x(x);

% строим график функции (по сеточным узлам)
figure(1)
plot(x,y,'-bo')
grid on
title('График функции y = f(x)','FontSize',TextSize)
xlabel('координата x','FontSize',TextSize)
ylabel('координата y','FontSize',TextSize)

% вводим расчетную точку и проверяем попадает ли она в диапазон
x0 = x_left - 1;
while  x0 > x_right || x0 < x_left
    x0_str =  inputdlg(['Введите значение х в даипазоне от:' num2str(x_left) ' до ' num2str(x_right)], 'Расчетная точка', 1, {num2str((x_right + x_left)/2*0.56)});
    x0 = str2double(x0_str{1});
end

% находим узлы слева и справа
k = find( (x > x0), 1, 'first'); % узел справа от точки
m = find( ( x < x0), 1, 'last'); % узел слева от точки, 

% -------------------------------------------------------------------------
% Интерполяция по Лагранжу первого порядка (линейная интерполяция)
yL1 = y(k-1)*(x0-x(k))/(x(k-1)-x(k)) + y(k)*(x0-x(k-1)) / (x(k)-x(k-1));
% погрешность
epsL1 = abs((yL1-f_ot_x(x0))/f_ot_x(x0))*100;

% Интерполяция по Лагранжу второго порядка (слева)
 yL2l = y(m - 1) * (x0 - x(m)) * (x0 - x(m + 1)) / (x(m - 1) - x(m)) / (x(m  - 1) - x(m + 1)) +...
     y(m) * (x0 - x(m - 1)) * (x0 - x(m + 1)) / (x(m) - x(m - 1)) / (x(m) - x(m + 1)) +...
     y(m + 1) * (x0 - x(m - 1)) * (x0 - x(m)) / (x(m + 1) - x(m - 1)) / (x(m + 1) - x(m));
 epsL2l = abs((yL2l - f_ot_x(x0)) / f_ot_x(x0))*100;

% Интерполяция по Лагранжу второго порядка (справа)
 yL2r = y(k - 2) * (x0 - x(k - 1)) * (x0 - x(k)) / (x(k - 2) - x(k - 1)) / (x(k - 2) - x(k)) +...
     y(k - 1) * (x0 - x(k - 2)) * (x0 - x(k)) / (x(k - 1) - x(k - 2)) / (x(k - 1) - x(k)) +...
     y(k) * (x0 - x(k - 2)) * (x0 - x(k - 1)) / (x(k) - x(k - 2)) / (x(k) - x(k - 1));

 epsL2r = abs((yL2r - f_ot_x(x0)) / f_ot_x(x0))*100;

% Интерполяция по Лагранжу третьего порядка
 yL3 = y(k - 3) * (x0 - x(k - 2)) * (x0 - x(k - 1)) * (x0 - x(k)) / (x(k - 3) - x(k - 2)) / (x(k - 3) - x(k - 1)) / (x(k - 3) - x(k)) +...
          y(k - 2) * (x0 - x(k - 3)) * (x0 - x(k - 1)) * (x0 - x(k)) / (x(k - 2) - x(k - 3)) / (x(k - 2) - x(k - 1)) / (x(k - 2) - x(k)) +...
          y(k - 1) * (x0 - x(k - 3)) * (x0 - x(k - 2)) * (x0 - x(k)) / (x(k - 1) - x(k - 3)) / (x(k - 1) - x(k - 2)) / (x(k - 1) - x(k)) +...
          y(k) * (x0 - x(k - 3)) * (x0 - x(k - 2)) * (x0 - x(k - 1)) / (x(k) - x(k - 3)) / (x(k) - x(k - 2)) / (x(k) - x(k - 1));

 epsL3 = abs((yL3 - f_ot_x(x0)) / f_ot_x(x0))*100;


% Сплайн-интерполяция
K = zeros(N - 1, 1); 
L =  zeros(N - 1, 1); 
K(1) = 0;
L(1) = 0;
for i=2:N-1
    hl = x(i) - x(i-1);
    hr = x(i+1) - x(i);
    a = hl/6;
    b = (hl+hr)/6;
    c = hr/6;
    d = (y(i+1)-y(i))/hr - (y(i)-y(i-1))/hl;
    K(i) = -c/(a*K(i-1)+2*b);
    L(i) = (d-a*L(i-1))/(a*K(i-1)+2*b);
end
m(N) = 0;
for i = (N-1):-1:1
    m(i) = K(i)*m(i+1)+L(i);
end
h = x(k) - x(k-1);
yS = m(k-1)*(x(k)-x0)^3*h/6 + m(k)*(x(k-1)-x0)^3*h/6 +...
    (y(k-1) - m(k-1)*h^2/6)*(x(k)-x0)/h + (y(k) - m(k)*h^2/6)*(x0-x(k-1))/h ;
epsS = abs((yS - f_ot_x(x0)) / f_ot_x(x0))*100;

% -------------------------------------------------------------------------
% вывод результатов 
hMsg = msgbox(  {'Результаты расчетов', ...
                ['Значение функции: в точке x0 = ' num2str(x0) ' y = ' num2str(f_ot_x(x0))], ...
                ['Интерполяция по Лагранжу 1 порядка y = ', num2str(yL1), ' погрешность eps = ', num2str(epsL1)], ...
                ['Интерполяция по Лагранжу 2 порядка справа y = ', num2str(yL2r), ' погрешность eps = ', num2str(epsL2r)], ...
                ['Интерполяция по Лагранжу 2 порядка слева y = ', num2str(yL2l), ' погрешность eps = ', num2str(epsL2l)], ...
                ['Интерполяция по Лагранжу 3 порядка y = ', num2str(yL3), ' погрешность eps = ', num2str(epsL3)], ...
                ['Cплайн - интерполяция y = ', num2str(yS),   ' погрешность eps = ', num2str(epsS)]
                }, ...
                'Результаты расчетов','warn');






%библиотечные функции

%Линейная интерполяция Matlab

y1 = interp1(x,y, x0, 'linear');
Eps1 = abs((y1 - f_ot_x(x0)) / f_ot_x(x0)) * 100;
yerror1 = abs((y1 - yL1) / yL1) * 100;
epserror1 = abs((Eps1 - epsL1) / epsL1) * 100;
fprintf('Сравнение линейной интерполяции \nОшибка по y: %d, по эпсилон: %d ', yerror1, epserror1);


%Интерполяция кубическим сплайном

y2 = interp1(x,y, x0, 'spline');
Eps2 = abs((y2 - f_ot_x(x0)) / f_ot_x(x0)) * 100;
yerror2 = abs((y2 - yS) / yS) * 100;
epserror2 = abs((Eps2 - epsS) / epsS) * 100;
fprintf('\n\nСравнение сплайн интерполяции \nОшибка по y: %d, по эпсилон: %d ', yerror2, epserror2);

%интерполяция по соседним точкам

y3 = interp1(x,y, x0, 'nearest');
Eps3 = abs((y3 - f_ot_x(x0)) / f_ot_x(x0)) * 100;
yerror3 = abs((y3 - yL2r) / yL2r) * 100;
epserror3 = abs((Eps3 - epsL2r) / epsL2r) * 100;
fprintf('\n\nСравнение интерполяции по Лагранжу второго порядка \nОшибка по y: %d, по эпсилон: %d ', yerror3, epserror3);

%интерполяция кубическим полиномом

y4 = interp1(x,y, x0, 'v5cubic');
Eps4 = abs((y4 - f_ot_x(x0)) / f_ot_x(x0)) * 100;
yerror4 = abs((y4 - yL3) / yL3) * 100;
epserror4 = abs((Eps4 - epsL3) / epsL3) * 100;
fprintf('\n\nСравнение интерполяции по Лагранжу третьего порядка \nОшибка по y: %d, по эпсилон: %d ', yerror4, epserror4);


%более подробная сетка

x_new =  x_left: h/2 : x_right;
yL1_new = zeros(1, length(x_new));
yL2r_new = zeros(1, length(x_new));
yL3_new = zeros(1, length(x_new));


% Интерполяция по Лагранжу первого порядка на подробной сетке

for i = 1:length(x_new)
    k_new = find( (x > x_new(i)), 1, 'first');
if (k_new - 1 < 1) || (k_new > length(x))

continue;
else
    yL1_new(i) = y(k_new - 1)*(x_new(i) - x(k_new)) / (x(k_new - 1) - x(k_new)) + y(k_new)*(x_new(i) - x(k_new - 1)) / (x(k_new) - x(k_new - 1));
end
end




% Интерполяция по Лагранжу второго порядка на подробной сетке

for i = 1:length(x_new)
    k_new = find( (x > x_new(i)), 1, 'first');
if (k_new - 1 < 1) || (k_new > length(x)) || (k_new - 2 < 1)

continue;
else
    yL2r_new(i) = y(k_new - 2) * (x_new(i) - x(k_new - 1)) * (x_new(i) - x(k_new)) / (x(k_new - 2) - x(k_new - 1)) / (x(k_new - 2) - x(k_new)) +...
     y(k_new - 1) * (x_new(i) - x(k_new - 2)) * (x_new(i) - x(k_new)) / (x(k_new - 1) - x(k_new - 2)) / (x(k_new - 1) - x(k_new)) +...
     y(k_new) * (x_new(i) - x(k_new - 2)) * (x_new(i) - x(k_new - 1)) / (x(k_new) - x(k_new - 2)) / (x(k_new) - x(k_new - 1));
end
end




% Интерполяция по Лагранжу третьего порядка на подробной сетке

for i = 1:length(x_new)
    k_new = find( (x > x_new(i)), 1, 'first');
if (k_new - 1 < 1) || (k_new > length(x)) || (k_new - 2 < 1) || (k_new - 3 < 1)

continue;
else
    yL3_new(i) = y(k_new - 3) * (x_new(i) - x(k_new - 2)) * (x_new(i) - x(k_new - 1)) * (x_new(i) - x(k_new)) / (x(k_new - 3) - x(k_new - 2)) / (x(k_new - 3) - x(k_new - 1)) / (x(k_new - 3) - x(k_new)) +...
          y(k_new - 2) * (x_new(i) - x(k_new - 3)) * (x_new(i) - x(k_new - 1)) * (x_new(i) - x(k_new)) / (x(k_new - 2) - x(k_new - 3)) / (x(k_new - 2) - x(k_new - 1)) / (x(k_new - 2) - x(k_new)) +...
          y(k_new - 1) * (x_new(i) - x(k_new - 3)) * (x_new(i) - x(k_new - 2)) * (x_new(i) - x(k_new)) / (x(k_new - 1) - x(k_new - 3)) / (x(k_new - 1) - x(k_new - 2)) / (x(k_new - 1) - x(k_new)) +...
          y(k_new) * (x_new(i) - x(k_new - 3)) * (x_new(i) - x(k_new - 2)) * (x_new(i) - x(k_new - 1)) / (x(k_new) - x(k_new - 3)) / (x(k_new) - x(k_new - 2)) / (x(k_new) - x(k_new - 1));
end
end


%Вывод графиков функции с графиком на грубой сетке

figure(2)
plot(x_new, yL1_new, '-bo', x, y, '-r.');
grid on
title('yL1(x)','FontSize', TextSize)
xlabel('координата x','FontSize',TextSize)
ylabel('координата y','FontSize',TextSize)



figure(3)
plot(x_new, yL2r_new, '-bo', x, y, '-r.');
grid on
title('yL2r(x)','FontSize', TextSize)
xlabel('координата x','FontSize',TextSize)
ylabel('координата y','FontSize',TextSize)


figure(4)
plot(x_new, yL3_new, '-bo', x, y, '-r.');
grid on
title('yL3(x)','FontSize', TextSize)
xlabel('координата x','FontSize',TextSize)
ylabel('координата y','FontSize',TextSize)



