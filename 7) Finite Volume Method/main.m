clc
close all
clear


TextSize = 20;
% Входные параметры
k = 1.4; % показатель адиабаты
Cp = 1006.43;
R = (k-1)/k*Cp;
p_inf = 1;
T_inf = 1;
M_inf = 3;
T0_inf = T_inf/tau_ot_M(M_inf,k);
V_inf = V_ot_M_T0_R_k(M_inf,T0_inf,R,k);


% load finite element mesh
pdeToolMesh = load('mesh.mat');
% convert and sort
[points,faces,leftCell,rightCell,numIntFace,numCell,numFaceInBc,startFaceInBc] = convertMeshElementBasedToFaceBsed(pdeToolMesh.t,pdeToolMesh.e,pdeToolMesh.p);
elementsByPoints = convertMeshFaceBsedToElementBased(faces,leftCell,rightCell,numIntFace,numCell);
numFace = length(leftCell);
numBcFace = numFace - numIntFace;


% расчет геометрических параметров сетки Vol,nx,ny, Sf, dx
[Sf,nx,ny,xf,yf,Vol,xC,yC,dx] = solveMeshData(points,faces,leftCell,rightCell,numIntFace,numCell);
return

%%
% создаем поля переменных + фиктивные ячейки + инициализация
[u1,u2,u3,u4] = initializationField(p_inf,V_inf,0.0,T_inf,k,R,numCell+numBcFace);
 %[u1,u2,u3,u4] = initializationField(p_inf,0.0,0.0,T_inf,k,R,numCell+numBcFace);

% добавляем индексы фиктивных ячеек как парвых 
rightCell((numIntFace+1):numFace) = ((numIntFace+1):numFace) - numIntFace + numCell;

% указываем список граничных условий (их типы), номер можно посмотреть в
% тулбоксе где генерировалась сетка 'SuperSonicInlet' 'SuperSonicOutlet' 'wall'
bcType = cellstr(char('SuperSonicInlet','wall',  'wall','wall','SuperSonicOutlet','wall'));
numBc = length(numFaceInBc);

% выводим сетку
figure(1)
plotUnstructuredDataAndMesh(points,faces',elementsByPoints,numFaceInBc,startFaceInBc,[],'mesh');
return
%%
% -------------------------------------------------------------------------

NumTimeStep =1000; % число шагов по времени
plot_interval=10; % интервал вывода графиков в шагах
CFL =0.75;
% 
t(1) = 0;
plot_time=plot_interval;
for i=1:NumTimeStep
    i
    % Добавляем граничные условия в фиктивных ячейках
    for bc = 1:numBc
        type = char(bcType(bc));
        switch (type)
            case {'SuperSonicOutlet'}
                [u1,u2,u3,u4] = addSuperSonicOutletBc(u1,u2,u3,u4,leftCell,rightCell,numFaceInBc(bc),startFaceInBc(bc));
            case {'SuperSonicInlet'}
                [u1,u2,u3,u4] = addSuperSonicInletBc(u1,u2,u3,u4,leftCell,rightCell,numFaceInBc(bc),startFaceInBc(bc),p_inf,V_inf,0.0,T_inf,R,k);
            case {'wall'}
                [u1,u2,u3,u4] = addWallBc(u1,u2,u3,u4,leftCell,rightCell,numFaceInBc(bc),startFaceInBc(bc));
            otherwise
                [u1,u2,u3,u4] = addWallBc(u1,u2,u3,u4,leftCell,rightCell,numFaceInBc(bc),startFaceInBc(bc));
        end
    end
    
    % расчитываем временной шаг 
    rhoV_mod = (u2.^2 + u3.^2).^0.5;
    dt= SolveTimeStepExplicitEiler(u1(1:numCell),rhoV_mod(1:numCell),u4(1:numCell),dx/2,CFL,k);
    t(i+1) = t(i) + dt;
    
    % расчитываем потоки
    % поворот системы координат
    [v1L,v2L] = rotateVectorToNormal(u2(leftCell),u3(leftCell),nx,ny);
    [v1R,v2R] = rotateVectorToNormal(u2(rightCell),u3(rightCell),nx,ny);
    % функция расчета потоков
  %[F1,F2,F3,F4,iteration] = SolveFluxExacRieman2D(u1(leftCell),v1L,v2L,u4(leftCell),...
                               % u1(rightCell),v1R,v2R,u4(rightCell),k);
     [F1,F2,F3,F4] = SolveFluxRoe2D(u1(leftCell),v1L,v2L,u4(leftCell),...
                                u1(rightCell),v1R,v2R,u4(rightCell),k);

    F1 = F1'; 
    F2 = F2';
    F3 = F3';
    F4 = F4';
    % поворот обратно
    [F2,F3] = rotateVectorToCartesian(F2,F3,nx,ny);
    % выполняем шаг по времени
    [u1,u2,u3,u4]=solveEvolution2D(F1.*Sf,F2.*Sf,F3.*Sf,F4.*Sf,u1,u2,u3,u4,leftCell,rightCell,dt,Vol,numIntFace);

 if i==plot_time
    figure(2)
    rho = u1(1:numCell);
    vx = u2(1:numCell)./rho;
    vy = u3(1:numCell)./rho;
    p = (k-1)*(u4(1:numCell) - 0.5*(vx.^2+vy.^2).*rho);
    T = p./(R*rho);
    %plotUnstructuredDataAndMesh(points,faces',elementsByPoints,numFaceInBc,startFaceInBc,p,'data');
     plotUnstructuredDataAndMesh(points,faces',elementsByPoints,numFaceInBc,startFaceInBc,rho,'data');
   %plotUnstructuredDataAndMesh(points,faces',elementsByPoints,numFaceInBc,startFaceInBc,(vx.^2+vy.^2).^0.5,'data');
   %plotUnstructuredDataAndMesh(points,faces',elementsByPoints,numFaceInBc,startFaceInBc,T,'data');
    plot_time=plot_time+plot_interval;
    pause
 end

end
    




















