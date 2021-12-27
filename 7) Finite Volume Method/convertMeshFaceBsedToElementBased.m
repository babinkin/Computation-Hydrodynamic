% Функция конвертирует из face-based неструктурированной сетке в
% element-based сетку
%
%
%
%
% Входные параметры:
%
% faces                 массив описывающий грани 
% leftCell              массив содержит левые ячейки к заданной грини
% rightCell             правые ячейки
% numIntFace            число внутренних граней
% numCell               количество внутренних ячеек
%
%
%
% Выходные параметры:
% 
% elementsByPoints      массив описывает элемент набором точек (узлов)
%
% =========================================================================
function elementsByPoints = convertMeshFaceBsedToElementBased(faces,leftCell,rightCell,numIntFace,numCell)
% --------------------------
numFace = length(leftCell);
numFaceInCell = 3; % FIXME ячейка может быть и не только треугольником
elementsByFace = zeros(numCell,numFaceInCell);
for f=1:numFace
    L = leftCell(f);
    R = rightCell(f);
    ind = find( elementsByFace(L,:)==0, 1,'first' );
    elementsByFace(L,ind) = f;
    if(f<=numIntFace)
        ind = find( elementsByFace(R,:)==0, 1,'first' );
        elementsByFace(R,ind) = -f;
    end
end
% --------------------------
elementsByPoints = zeros(numCell,3); 
for c=1:numCell
    numFaceInCell = 3;
    for f=1:(numFaceInCell-1)
        face = elementsByFace(c,f);
        if (face>0) % left face
            p1 = faces(face,1);
            p2 = faces(face,2);
        else
            p1 = faces(-face,2);
            p2 = faces(-face,1);
        end
        % ---------
        indP = find( elementsByPoints(c,:)==p1, 1,'first' );
        if(isempty(indP))
            ind0 = find( elementsByPoints(c,:)==0, 1,'first' );
            elementsByPoints(c,ind0) = p1;
        end
        indP = find( elementsByPoints(c,:)==p2, 1,'first' );
        if(isempty(indP))
            ind0 = find( elementsByPoints(c,:)==0, 1,'first' );
            elementsByPoints(c,ind0) = p2;
        end
    end
end


end