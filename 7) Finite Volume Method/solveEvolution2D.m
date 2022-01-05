% функция выполняет шаг интегрирования по времени
%
%
%
%
% =========================================================================
function [u1,u2,u3,u4]=solveEvolution2D(F1,F2,F3,F4,u1,u2,u3,u4,leftCell,rightCell,dt,Vol,numIntFace)

numFace = length(F1);
for f=1:numFace
    left = leftCell(f);
    right = rightCell(f);
    %
    tL = dt/Vol(left);
    u1(left) = u1(left) - tL*F1(f);
    u2(left) = u2(left) - tL*F2(f);
    u3(left) = u3(left) - tL*F3(f);
    u4(left) = u4(left) - tL*F4(f);
    %
    if (f<=numIntFace)
        tR = dt/Vol(right);
        u1(right) = u1(right) + tR*F1(f);
        u2(right) = u2(right) + tR*F2(f);
        u3(right) = u3(right) + tR*F3(f);
        u4(right) = u4(right) + tR*F4(f);
    end
    
end



end