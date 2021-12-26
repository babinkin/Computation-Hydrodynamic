% Функция конвертирует из element-based неструктурированной сетке в
%  face-based сетку
%
%
%
%
% Выходные параметры:
%
% points                массив координат точек
% faces                 массив описывающий грани 
% leftCell              массив содержит левые ячейки к заданной грини
% rightCell             правые ячейки
% numIntFace            число внутренних граней
% numCell               количество внутренних ячеек
% numFaceInBc           количество граней в границе
% startFaceInBc         стартовый номер граничной гране в общем списке
%
%
%
% Входные параметры: импортируют данные о сетке из pdetool
% 
% elementsByPoints      массив описывает элемент набором точек (узлов)
% boundary              массив описывающий границу
% p                     массив точеу
%
% =========================================================================
function [points,faces,leftCell,rightCell,numIntFace,numCell,numFaceInBc,startFaceInBc] = convertMeshElementBasedToFaceBsed(elementsByPoints,boundary,p)

% -------------------------------------------------------------------------
% convert to face-based notation
% 
% for all triang
numCell = length(elementsByPoints(1,:));
faces = [];
leftCell = [];
rightCell = [];
for tri=1:numCell
    tri;
    numFaceInCell = 3;
    for f=1:(numFaceInCell-1)
        p1 = elementsByPoints(f,tri);
        p2 = elementsByPoints(f+1,tri);
        % find right (ищем есть ли в массиве faces грани с такими же узлами
        % но в обр порядке, может быть только не более 1 для 2Д сеток, и
        % ячейка с такой гранью будет соседнй к данной)
        if(length(faces)>0)
            indFaceRightCell = find( faces(:,1)==p2 & faces(:,2)==p1 );
        else
            indFaceRightCell=[];
        end
        indFaceRightCell;
        % если не нашли, то зиписываем и устанавливаем  leftCell
        % если нашли, то значение в массиве leftCell для этой грани уже
        % установленно верно и устанавливаем для этой грани rightCell как
        % tri. Эта грань будет внутренней.
        if(indFaceRightCell)
            rightCell(indFaceRightCell) = tri;
        else
            % добавляем грать p1-p2 в массив faces, в массив левых ячеек
            % добавляем индекс ячейки i
            faces(end+1,1) = p1; % push_back
            faces(end,2) = p2; % set
            leftCell(end+1) = tri; % push_back
        end
    end
    % то же для послед грани в ячейке
    p1 = elementsByPoints(numFaceInCell,tri);
    p2 = elementsByPoints(1,tri);
    indFaceRightCell = find( faces(:,1)==p2 & faces(:,2)==p1 );
    if(indFaceRightCell)
        rightCell(indFaceRightCell) = tri;
    else
        faces(end+1,1) = p1; % push_back
        faces(end,2) = p2; % set
        leftCell(end+1) = tri; % push_back
    end
end
% там где rightCell Осталось == 0 это граничные грани, размер rightCell
% делаем такой же как и у leftCell
if(length(rightCell)<length(leftCell))
    rightCell(length(leftCell)) =0;
end

% -------------------------------------------------------------------------
% сортируем грани сначала внутренние, потом граничные 1 Гу, потом 2 го и
% т.д. !!!! FIXME сделать сортировку с алгоритмом уменьшения ширины ленты
% ----------
% подсчитываем сколько всего ГУ и сколько граней в каждом
numBc = max(boundary(5,:));
numFace = length(faces(:,1));
numBcFace = length(boundary(5,:));
numIntFace = numFace - numBcFace;
startFaceInBc(1) = numIntFace + 1;
for bc = 1:numBc-1
    ind  = find(boundary(5,:)==bc);
    size = length(ind);
    numFaceInBc(bc) = size;
    startFaceInBc(bc+1) = startFaceInBc(bc) + numFaceInBc(bc);
end
size  = length(find(boundary(5,:)==numBc));
numFaceInBc(numBc) = size;

% Sort face
for bc = 1:numBc
    ind  = find(boundary(5,:)==bc);
    for f = 1:numFaceInBc(bc)
        indBc = ind(f);
        if(boundary(7,indBc) == 1)
            p1 = boundary(2,indBc);
            p2 = boundary(1,indBc);
        else
            p1 = boundary(1,indBc);
            p2 = boundary(2,indBc);
        end
        
        
        indFace1 = find( faces(:,1)==p1 & faces(:,2)==p2 );
%         rightCell(indFace1)
        indFace2 = startFaceInBc(bc) + f - 1;
        % обмен
        L = leftCell(indFace1);
        R = rightCell(indFace1);
        faces(indFace1,1) = faces(indFace2,1);
        faces(indFace1,2) = faces(indFace2,2);
        leftCell(indFace1) = leftCell(indFace2);
        rightCell(indFace1) = rightCell(indFace2);
        %
        faces(indFace2,1) = p1;
        faces(indFace2,2) = p2;
        leftCell(indFace2) = L;
        rightCell(indFace2) = R;
    end
end

% здесь нужно бы реализовать сортировку ячеек и граней для уменьшения
% ширины ленты и "порядок верхнего треугольника" up-triang order
%


% -----------------------------
points = p';
leftCell = leftCell';
rightCell = rightCell';
numFaceInBc = numFaceInBc';
startFaceInBc = startFaceInBc';


end