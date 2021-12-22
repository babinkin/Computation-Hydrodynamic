% функция заносит в фиктивные ячейки параметры соответсвующие граничному
% условию
%
%
%
% 
% =========================================================================

function [u1,u2,u3,u4] = addWallBc(u1,u2,u3,u4,leftCell,rightCell,num,start)
for f=start:(start+num -1)
    left = leftCell(f);
    right = rightCell(f);
    u1(right) = u1(left);
    u2(right) = -u2(left);
    u3(right) = -u3(left);
    u4(right) = u4(left);
end

end











