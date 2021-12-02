clc
close
clear

TextSize = 20;
%==========================================================================
% Количество узлов 
N = 100;
% показатель адиабаты
k = 1.4;

var = 5;
% координаты точек, необходимые для построения контура сопла
x0 = 0;
x_zv = 1;
x_a = 3;
r0 = 1.5;
r_zv = 1;
r_a = 2;

p0 = (10 + 3 * var) * 101325;
T0 = 300;

% строим контур
x = x0:(x_a - x0)/(N-1):x_a;
% не эффективный алгоритм зато векторный
cond = x<x_zv;
r =  (r0 + (r_zv-r0)/(x_zv -x0)*(x-x0) ).*cond + (r_zv + (r_a-r_zv)/(x_a-x_zv)*(x - x_zv)).*(1-cond);

% контур сопла
figure(1)
plot(x,r,'-bo')
grid on
title(['профиль сопла'],'FontSize',TextSize)
xlabel('координата x','FontSize',TextSize)
ylabel('радиус сечения r','FontSize',TextSize)

F_zv = pi*r_zv^2; % площадь критического сечения
F = pi*r.^2; % площади всех сечений
q = F_zv./F; % функция приведенного расхода

% решение трансцендентного уравнения
% нахождение числа Маха

M = zeros(1, N);
Pi = zeros(1, N);
Eps = zeros(1, N);
Tau = zeros(1, N);
 for i=1:N
     if x(i)<1
         M_start = 0.1; % Начальное приближение в дозвуковой части
     else
         M_start = 1.5; % Начальное приближение в сверхвуковой части
     end
     M(i) = SolveQFun(q(i),k,M_start); 
     % по формулам изоэнтропического течения
    Pi(i) = (1 + ((k - 1) / 2) * M(i)^2)^(-k / (k - 1));
    Eps(i) = (1 + ((k - 1) / 2) * M(i)^2)^(-1 / (k - 1));
     Tau(i) = (1 + ((k - 1) / 2) * M(i)^2)^(-1);
 end

% строим графики результатов 

figure(2);
plot(x, M, '-ro');
title('Графики зависимостей', 'FontSize', TextSize);
xlabel('x', 'FontSize', TextSize);
ylabel('Параметр', 'FontSize', TextSize);
grid on;
hold on;
plot(x, Pi, '-bo');
plot(x, Eps, '-go');
plot(x, Tau, '-yo');
plot(x, q, '-co');
hold off;
legend({'Mаха', 'давлений', 'плотностей', 'температур', 'расхода'}, 'Location', 'best');


Mol = 0.029; % молярная масса воздуха (кг / моль)
R_gaz = 8.314; % универсальная газовая постоянная (Н * м / моль * К)
ro0 = 15.341; % плотность торможения воздуха (кг / м^3)
v_max = sqrt(2 * T0 * (k / (k - 1)) * R_gaz / Mol);
p = p0.*Pi; % давление
rho = ro0.*Eps; % плотность
T = T0.*Tau; % температура
v = zeros(1, N); % скорость
for i = 1:N
    v(i) = v_max * (((k - 1) / 2) * M(i)^2 / (1 + (((k - 1) / 2) * M(i)^2)))^(1 / 2);
end


figure(3);
plot(x, p, '-ro');
title('Распределение давления', 'FontSize', TextSize);
xlabel('x', 'FontSize', TextSize);
ylabel('Па', 'FontSize', TextSize);
grid on;

figure(4);
plot(x, rho, '-ro');
title('Распределение плотности', 'FontSize', TextSize);
xlabel('x', 'FontSize', TextSize);
ylabel('кг / м^3', 'FontSize', TextSize);
grid on;

figure(5);
plot(x, T, '-ro');
title('Распределение температуры', 'FontSize', TextSize);
xlabel('x', 'FontSize', TextSize);
ylabel('К', 'FontSize', TextSize);
grid on;

figure(6);
plot(x, v, '-ro');
title('Распределение скорости', 'FontSize', TextSize);
xlabel('x', 'FontSize', TextSize);
ylabel('м / с', 'FontSize', TextSize);
grid on;


figure(7);
plot(x, q, '-ro');
title('Распределение расхода', 'FontSize', TextSize);
xlabel('x', 'FontSize', TextSize);

grid on;


p_atm = 101325; % нормальное атмосферное давление (Па)
n = p(N) / p_atm;
fprintf('Число нерасчетности истечения:%f. ', n);



M_a = M(N);
fprintf('Геометрическое число Маха сопла: %f.', M_a);











