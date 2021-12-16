clc
close
clear

TextSize = 20;

Gamma = 1.4; % показатель адиабаты
Cp = 1006.43;
R = (Gamma - 1) / Gamma * Cp;

test = 1; % номер теста
[rhoL, uL, pL, rhoR, uR, pR, x_dis, NumCell, t_fin, CFL, ...
    NumTimeStep, BC_L, BC_R] = ShockTubeInitData(test);

TL = pL / (R * rhoL);
TR = pR / (R * rhoR);

plot_interval = 1; % интервал вывода графиков в шагах
NumSmallStep = 0; % количество шагов разгона

OrderInSpace = 1; % порядок аппроксимации по пространству
RiemannSolver = 'Roe'; % 'Exac', 'Roe'
PlotType = 'Rho'; % 'Rho', 'Vel', 'Press', 'IntEnerg', 'Temp', 'RhoV'

% Сетка

nu = 1 ; % коэффициент симметри
x_left = 0; % левая граница
x_right = 1; % правая граница
delta_h = (x_right - x_left) / NumCell;
x_i = (x_left+0.5*delta_h):delta_h:x_right; % центры ячеек

ind_x_dis = floor(x_dis / delta_h); % количество ячеек, поместившихся целиком в x_shock

cut = zeros(1, NumCell);

for i = 1:NumCell
   
    if x_i(i) < 0.335
        
        cut(i) = 1;
   
    elseif x_i(i) > 0.675

        cut(i) = 0;

    else
        
    end
    
end

%Переменные.

U1old = zeros(1, NumCell+2);
U2old = zeros(1, NumCell+2);
U3old = zeros(1, NumCell+2);

pLeftInit = 2:ind_x_dis+1;
pRightInit = ind_x_dis+2:NumCell+1;

U1old(pLeftInit) = rhoL;
U1old(pRightInit) = rhoR;

U2old(pLeftInit) = rhoL * uL;
U2old(pRightInit) = rhoR * uR;

U3old(pLeftInit) = pL / (Gamma - 1) + 0.5 * uL^2 * rhoL;
U3old(pRightInit) = pR / (Gamma - 1) + 0.5 * uR^2 * rhoR;

% Граничные условия:
U1old(1) = U1old(2);
U2old(1) = BC_L * U2old(2);
U3old(1) = U3old(2);
U1old(end) = U1old(end-1);
U2old(end) = BC_R * U2old(end-1);
U3old(end) = U3old(end-1);

U1new = U1old;
U2new = U2old;
U3new = U3old;

% Временные переменные.
t(1) = 0;
plot_time = plot_interval;
% Цикл тестов.

for i = 1:NumTimeStep
    
    dt = SolveTimeStepExplicitEiler(U1old, U2old, U3old, delta_h, CFL, Gamma);
    
if i <= NumSmallStep
    
    dt = i / NumSmallStep * dt;
    
end

t(i+1) = t(i) + dt;

switch (OrderInSpace)
    
    case {1}
        
        switch (RiemannSolver)
            
            case {'Exac'}
                
                [F1, F2, F3, iterationExac] = SolveFluxExacRieman(U1old(1:end-1), U2old(1:end-1), U3old(1:end-1), U1old(2:end), U2old(2:end), U3old(2:end), Gamma); % Godunov
              
                if iterationExac > 1
                    
                    isp(iterationExac);

                end
                
            case {'Roe'}
               
                [F1, F2, F3] = SolveFluxRoe(U1old(1:end-1), U2old(1:end-1), U3old(1:end-1), U1old(2:end), U2old(2:end), U3old(2:end), Gamma); % Roe
            
            otherwise
                
                [F1, F2, F3, iterationExac] = SolveFluxExacRieman(U1old(1:end-1), U2old(1:end-1), U3old(1:end-1), U1old(2:end), U2old(2:end), U3old(2:end), Gamma); % Godunov

        end
        
end

if OrderInSpace > 0
  
    [U1new, U2new, U3new] = SolveEvolutionExplFirstOrder(F1, F2, F3, U1old, U2old, U3old, nu, dt, delta_h, 0, Gamma);

end

U1old = U1new;
U2old = U2new;
U3old = U3new;

% Рассчитываем граничные условия:
U1old(1) = U1old(2); U1old(end) = U1old(end-1);
U2old(1) = BC_L * U2old(2); U2old(end) = BC_R * U2old(end-1);
U3old(1) = U3old(2); U3old(end) = U3old(end-1);

% Конец цикла
if i == plot_time
       
    PlotTestGraphics(1, U1new, U2new, U3new, rhoL, uL, pL, rhoR, uR, pR, x_i, t, x_dis, Gamma, R, TextSize, test, PlotType);
        plot_time = plot_time + plot_interval;
        
end

% Выход, если досчитали:
    if t(i+1) > t_fin
       
        PlotTestGraphics(1, U1new, U2new, U3new, rhoL, uL, pL, rhoR, uR, pR, x_i, t, x_dis, Gamma, R, TextSize, test, PlotType);
        disp(['Количество итераций: ', num2str(i), '.']);
       
        break;

    end
    
end











































