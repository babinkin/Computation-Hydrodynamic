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
%
k = 5; % номер варианта
n = k * 5 + 100; % количество неизвестных
% создаем трехдиаганальную матрицу 
A = gallery('lesp', n );
% A = full(gallery('tridiag',n,-1,4,-1));
b_right = k * ones(1, n)' ;% вектор правых частей
% Находи L1 норму
L1 =  max(sum(abs(A)));
% L1 = norm(A,1); % библиотечная функция
cond1 = L1*norm(inv(A),1);
% L_inf норма
L_inf = max(sum(abs(A')));
% число обусловленности
cond_inf = L_inf * norm(inv(A), 1);

fprintf('L1 норма: %f.\nЧисло обусловленности: %f.', L1 , cond1);
fprintf('\nLinf норма: %f.\nЧисло обусловленности: %f.', L_inf , cond_inf);


% -------------------------------------------------------------------------
% алгоритм прогонки
% прямой ход
% прогоночные коэффициенты

a = [0; diag(A, -1)];
b = diag(A);
c = [diag(A, 1); 0]; % массивы элементов трех диагоналей
alpha = zeros(n-1, 1, "like", b_right);
beta = zeros(n-1, 1, "like", b_right);
alpha(1) = c(1) / b(1);
beta(1) = b_right(1) / c(1); % инициализация прямого хода прогонки
for i = 2:(n-1)
    alpha(i) = c(i) / (b(i) - a(i) * alpha(i-1));
    beta(i) = (b_right(i) - a(i) * beta(i-1)) / (b(i) - a(i) * alpha(i-1));
end % формулы прямого хода
x = zeros(n, 1, "like", b_right);
x(n) = (b_right(n) - a(n) * beta(n-1)) / (b(n) - a(n) * alpha(n-1)); % инициализация обрат
for i = (n-1):-1:1
    x(i) = beta(i) - alpha(i) * x(i+1);
end % формулы обратного хода

fprintf('\nРешение методом подгонки');
fprintf('\n%f', x);

%метод Якоби
x_jac_old = ones(n, 1, "like", b_right);
x_jac_new = zeros(n, 1, "like", b_right); % начальное приближение
counter = 0; % подсчет количества итераций
while norm(x_jac_new - x_jac_old, 2) > 1e-6 % критерий сходимости
    x_jac_old = x_jac_new;
    for i = 1:n
        sigma = 0;
        for j = 1:n
            if j ~= i
                sigma = sigma + A(i,j) * x_jac_old(j);
            end
        end
        x_jac_new(i) = (1 / A(i,i)) * (b_right(i) - sigma);
   
    end
    counter = counter + 1;
end


fprintf('\nКоличество итераций метода Якоби: %d.', counter);


fprintf('\nРешение методом Якоби');
fprintf('\n%f', x_jac_new);
% 
% figure(1)
% imagesc(A);
% title('Структура матрицы', 'FontSize', TextSize)
% axis square;
% colorbar;



e = eig(A);
fprintf('\nСобственные числа: ');
fprintf('\n %d.', e);

L_2 = max(sqrt(eig(transpose(A) * A))); % L_2-норма
cond_2 = L_2 * max(sqrt(eig(transpose(inv(A)) * inv(A)))); % число обусловленности для L_2
L_E = norm(A, 2); % евклидова норма
cond_E = L_E * norm(inv(A), 2); % число обусловленности для евклидовой нормы

fprintf('\nL2 норма: %f.\nЧисло обусловленности: %f.', L_2 , cond_2);
fprintf('\nL_E норма: %f.\nЧисло обусловленности: %f.', L_E , cond_E);


%метод Гаусса матлаб
gauss = A \ b_right; 

%метод обратной матрицы
invmatrix = inv(A) *  b_right;

%LU разложение
[L, U] = lu(A);

y = L\b_right;
ludecomposition = U \ y;

%qr разложение

qrdecomposition = linsolve(A, b_right);


err_tridiag_gauss = gauss - x;
err_tridiag_invmatrix = invmatrix - x;
err_tridiag_ludecomposition = ludecomposition - x;
err_tridiag_qrdecomposition= qrdecomposition - x;


err_tridiag2_gauss = gauss - x_jac_new;
err_tridiag2_invmatrix = invmatrix - x_jac_new;
err_tridiag2_ludecomposition = ludecomposition - x_jac_new;
err_tridiag2_qrdecomposition= qrdecomposition - x_jac_new;



figure(1);
plot(x - x_jac_new, '-ro');
grid on;
title('Засисимость невязки от итераций', 'FontSize', TextSize);
ylabel('Невязка', 'FontSize', TextSize);
xlim([0 20]);

% figure(2);
% plot(err_tridiag_gauss, '  -r.');
% title('Cравнение ошибки метода Гаусса (матлаб))', 'FontSize', TextSize);
% ylabel('Ошибка', 'FontSize', TextSize);
% xlim([0 125]);
% ylim([-0.5 0.5]);
% hold on;
% plot(err_tridiag2_gauss, 'bo');
% hold off;
% legend({'C методом прогонки', 'С методом Якоби'}, 'Location', 'best');
% 
% figure(3);
% plot(err_tridiag_invmatrix, '  -r.');
% title('Cравнение ошибки метода обратной матрицы (матлаб))', 'FontSize', TextSize);
% ylabel('Ошибка', 'FontSize', TextSize);
% xlim([0 125]);
% ylim([-0.5 0.5]);
% hold on;
% plot(err_tridiag2_invmatrix, 'bo');
% hold off;
% legend({'C методом прогонки', 'С методом Якоби'}, 'Location', 'best');
% 
% figure(4);
% plot(err_tridiag_ludecomposition, ' -r.');
% title('Cравнение ошибки LU разложения  (матлаб) )', 'FontSize', TextSize);
% ylabel('Ошибка', 'FontSize', TextSize);
% xlim([0 125]);
% ylim([-0.5 0.5]);
% hold on;
% plot(err_tridiag2_ludecomposition, 'bo');
% hold off;
% legend({'C методом прогонки', 'С методом Якоби'}, 'Location', 'best');
% 
% figure(5);
% plot(err_tridiag_qrdecomposition, ' -r.');
% title('Cравнение ошибки QR разложения (матлаб))', 'FontSize', TextSize);
% ylabel('Ошибка', 'FontSize', TextSize);
% xlim([0 125]);
% ylim([-0.5 0.5]);
% hold on;
% plot(err_tridiag2_qrdecomposition, 'bo');
% hold off;
% legend({'C методом прогонки', 'С методом Якоби'}, 'Location', 'best');
% 


r_tridiag_gauss = b_right - A * gauss;
r_tridiag_invmatrix = b_right - A * invmatrix;
r_tridiag_ludecomposition = b_right - A * ludecomposition;
r_tridiag_qrdecomposition = b_right - A * qrdecomposition;
r_tridiag = b_right - A * x;



% figure(6);
% plot(r_tridiag, 'r.');
% grid on;
% title('График невязки метода прогонки', 'FontSize', TextSize);
% ylabel('Невязка', 'FontSize', TextSize);
% xlim([0 107]);
% ylim([-0.5 0.5]);
% 
% figure(7);
% plot(r_tridiag_gauss, 'r.');
% grid on;
% title('График невязки метода метода Гаусса (матлаб', 'FontSize', TextSize);
% ylabel('Невязка', 'FontSize', TextSize);
% xlim([0 107]);
% ylim([-0.5 0.5]);
% 
% figure(8);
% plot(r_tridiag_invmatrix, 'r.');
% grid on;
% title('График невязки метода обратной матрицы (матлаб)', 'FontSize', TextSize);
% ylabel('Невязка', 'FontSize', TextSize);
% xlim([0 107]);
% ylim([-0.5 0.5]);
% 
% figure(9);
% plot(r_tridiag_ludecomposition, 'r.');
% grid on;
% title('График невязки LU разложения  (матлаб)', 'FontSize', TextSize);
% ylabel('Невязка', 'FontSize', TextSize);
% xlim([0 107]);
% ylim([-0.5 0.5]);
% 
% figure(10);
% plot(r_tridiag_qrdecomposition, 'r.');
% grid on;
% title('График невязки QR разложения (матлаб)', 'FontSize', TextSize);
% ylabel('Невязка', 'FontSize', TextSize);
% xlim([0 107]);
% ylim([-0.5 0.5]);