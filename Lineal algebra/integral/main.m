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
%==========================================================================

% Количество узлов 
N = 100;
% Диапазон изменения x
x_left = 0;
x_right = 5;
x0 = 1.4;
% вычисляем шаг по x (равномерный)
h = (x_right - x_left) / (N - 1);

% получаем значение функции в узловых точках
x = x_left: h :x_right;
y = f_int_ot_x(x);

% строим график функции (по сеточным узлам)
figure(1)
plot(x,y,'-bo')
grid on
title( [ 'График функции y=f(x) ' ], 'FontSize', TextSize)
xlabel('координата x', 'FontSize', TextSize)
ylabel('координата y', 'FontSize', TextSize)

% -------------------------------------------------------------------------
% метод левых прямоугольников
I_left = sum(f_int_ot_x(x(1:N)))*h;
eps_left = max(d1_f(x))*h/2;

% метод центральных прямоугольников
x_c = x(1:N-1)+0.5*h;
I_cent = sum(f_int_ot_x(x_c))*h;
eps_cent = max(d2_f(x))*h^2/24;

% метод трапеций
 I_trap = h/2*(f_int_ot_x(x_left) + f_int_ot_x(x_right) + 2*sum(f_int_ot_x(x(2:N-1))));
 eps_trap = max(d2_f(x))*h^2/12;

% метод Симпсона
I_Simp = ( f_int_ot_x(x_left)+f_int_ot_x(x_right) + 4*sum(f_int_ot_x(x_c)) +2*sum(f_int_ot_x(x(2:N-1))) )*h/6;
 eps_Simp = max(d4_f(x))*h^4/2880;


% вычисляем производную разностью вправо
k=50;
dydx_1_right = ( y(k+1) - y(k) )/( x(k+1) - x(k) );
dydx_exac = d1_f(x(k));
eps_dydx = abs((dydx_1_right-dydx_exac)/dydx_exac)*100;



%произовдная со вторым порядком аппроксимации

dydx_2 = (f_int_ot_x(x0 + h) - f_int_ot_x(x0 - h)) / (h * 2);
eps_dydx2 = abs((dydx_2 - d1_f(x0)) / d1_f(x0)) * 100;

%интеграл при помощи библиотечных функций

syms x_matlab;
f  = cos(10.*x_matlab).*exp(-x_matlab/2);
int_matlab = vpaintegral(f, x_matlab, [x_left, x_right]);
eps_matlab = abs((int_matlab - I_cent) / I_cent) * 100;




% вывод результатов
disp(['метод левых прям.   ' num2str(I_left) '   погрешность   ' num2str(eps_left)])
disp(['метод центр прям    ' num2str(I_cent) '   погрешность   ' num2str(eps_cent)])
 disp(['метод трапеций      ' num2str(I_trap) '   погрешность   ' num2str(eps_trap)])
 disp(['метод Симпсона      ' num2str(I_Simp) '   погрешность   ' num2str(eps_Simp)])
disp(['Производная в в узле к = ' num2str(k) ' с первым порядком при  dydx= ' num2str(dydx_1_right) ' ошибка eps = ' num2str(eps_dydx)])
disp(['Производная со вторым порядком аппроксимации', num2str(dydx_2), ' погрешность ' ,  num2str(eps_dydx2)])
disp(['Интеграл с использованием матлаба',  num2str(double(int_matlab)) ' погрешность в сравнении с методом прямоугольников'  num2str(double(eps_matlab))])




















