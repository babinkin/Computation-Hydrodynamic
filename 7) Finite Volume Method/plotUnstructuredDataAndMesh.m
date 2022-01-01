% Функция отображает поле параметра на неструктурированной сетке
%
%
%
%
% Входные параметры:
%
% points                массив координат точек
% faces                 массив описывающий грани 
% leftCell              массив содержит левые ячейки к заданной грини
% rightCell             правые ячейки
% numIntFace            число внутренних граней
% numCell               количество внутренних ячеек

% elementsByPoints      сетка в element-based нотации
% numFaceInBc           количество граней в границе
% startFaceInBc         стартовый номер граничной гране в общем списке
% cdata                 массив параметров данных
% plotType              тип вывода 'mesh' - просто сетка 'data' данные
% 
%
% =========================================================================
function plotUnstructuredDataAndMesh(points,faces,elementsByPoints,numFaceInBc,startFaceInBc,cdata,plotType)
% -------------------------------------------------------------------------
% визуалиция сетки
switch (plotType)
    case {'data'}
        patch('Faces',elementsByPoints,'Vertices',points,'FaceVertexCData',cdata,'FaceColor','flat','EdgeColor','none','CDataMapping','scaled');
        colorbar
    case {'mesh'}
        patch('Faces',elementsByPoints,'Vertices',points,'FaceColor','w','EdgeColor','g');
    otherwise
        patch('Faces',elementsByPoints,'Vertices',points,'FaceColor','w','EdgeColor','g');
end
% ----------------------------------------
numBc = length(numFaceInBc);
for bc = 1:numBc
    indP1 = faces(1,startFaceInBc(bc):(startFaceInBc(bc)+numFaceInBc(bc) -1 ));
    x = mean(points(indP1,1));
    y = mean(points(indP1,2));
    text(x,y,['B' num2str(bc)],'FontSize',16)
end


axis equal tight


end