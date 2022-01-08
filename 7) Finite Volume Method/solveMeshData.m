% Функция расчитывает геометрические параметры сетки
%
%
% Входные параметры
%
% points                массив координат точек
% faces                 массив описывающий грани 
% leftCell              массив содержит левые ячейки к заданной грини
% rightCell             правые ячейки
% numIntFace            число внутренних граней
% numCell               количество внутренних ячеек
%
% Выходные параметры
%
% Sf                    площадь (длина) граней
% nx                    x компонента вектора нормали грани
% ny                    y компонента 
% xf                    центр грани
% yf                    
% Vol                   объем (площадь в 2Д случае) ячейки
% xC                    центр ячейки
% yC                    
% dx                    диаметр вписанной окружности в ячейку
% 
% =========================================================================
function [Sf,nx,ny,xf,yf,Vol,xC,yC,dx] = solveMeshData(points,faces,leftCell,rightCell,numIntFace,numCell)
% calc face normals
p1 = faces(:,1);
p2 = faces(:,2);
x1 = points(p1,1);
y1 = points(p1,2);
x2 = points(p2,1);
y2 = points(p2,2);
Sf = ((x2 - x1).^2 + (y2 - y1).^2).^0.5;
nx = (y2 - y1)./Sf;
ny = -(x2 - x1)./Sf;
xf = 0.5*(x2 + x1);
yf = 0.5*(y2 + y1);


% calc face centr
xC = zeros(numCell,1);
yC = zeros(numCell,1);
numFaceInCell = 3;
numFace = length(leftCell);
for f=1:numFace
    xC(leftCell(f)) = xC(leftCell(f)) +  xf(f);
    yC(leftCell(f)) = yC(leftCell(f)) +  yf(f);
    if(f<=numIntFace)
        xC(rightCell(f)) = xC(rightCell(f))  + xf(f);
        yC(rightCell(f)) = yC(rightCell(f))  + yf(f);
    end
end
xC = xC./numFaceInCell;
yC = yC./numFaceInCell;


% check normals
numModifNormal = 0;
for f=1:numFace
    dr_x = xf(f) - xC(leftCell(f));
    dr_y = yf(f) - yC(leftCell(f));
    dot = dr_x*nx(f) + dr_y*ny(f);
    if(dot<0)
        numModifNormal = numModifNormal +1;
        temp = faces(f,1);
        faces(f,1) = faces(f,2);
        faces(f,2) = temp;
        nx(f) = -nx(f);
        ny(f) = -ny(f);
    end    
end

% calc volume
Vol = zeros(numCell,1);
for f=1:numFace
    temp = xf(f)*nx(f)*Sf(f);
    Vol(leftCell(f))= Vol(leftCell(f)) + temp;
    if (f<=numIntFace)
        Vol(rightCell(f))= Vol(rightCell(f)) - temp;
    end
end



% расчитываем радиус вписанной окружности FIXME для начала берем как
% диаметр круга равной площади
dx = (4/pi*Vol).^0.5;


    
end















