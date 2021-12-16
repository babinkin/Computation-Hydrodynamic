function [rhoL, uL, pL, rhoR, uR, pR, x_dis, NumCell, t_fin, CFL, MaxIter, BC_L, BC_R] = ShockTubeInitData(test)

switch (test)
    
case {1}
% Test 1 (Modified Sod):
        x_dis = 0.3; % Position of diaphragm 1
        NumCell = 100; % Number of computing cells
        t_fin = 0.2; % Output time
        rhoL = 1.0; % Initial density  on left section of tube
        uL = 0.75; % Initial velocity on left section of tube
        pL = 1.0; % Initial pressure on left section of tube
        rhoR = 0.125; % Initial density  on right section of tube
        uR = 0; % Initial velocity on right section of tube
        pR = 0.1; % Initial pressure on right section of tube
        CFL = 0.9; % Courant number coefficient
        MaxIter = 10000000; % Maximum number of time steps
        BC_L = 1; % Type of left boundary conditions 1  upstream, -1 wall
        BC_R = -1; % Type of right boundary conditions

end




end

% Входные параметры:
% test - номер задачи.
%
% Выходные параметры:
% rhoL, uL, pL - значение параметров слева,
% rhoR, uR, pR - значение параметров справа,
% t_fin - время завершения,
% x_dis - место расположения разрыва,
% NumCell - количество ячеек.