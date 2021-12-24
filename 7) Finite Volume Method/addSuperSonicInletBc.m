% функция заносит в фиктивные ячейки параметры соответсвующие граничному
% условию addSuperSonicInletBc
%
%
%
% 
% =========================================================================

function [u1,u2,u3,u4] = addSuperSonicInletBc(u1,u2,u3,u4,leftCell,rightCell,num,start,p_inf,vx_inf,vy_inf,T_inf,R,k)
rho = p_inf/(R*T_inf);
rhoV = rho*vx_inf;
rhoU = rho*vy_inf;
rhoE = p_inf/(k-1) + 0.5*rho*(vx_inf*vx_inf + vy_inf*vy_inf);

for f=start:(start+num -1)
    right = rightCell(f);
    u1(right) = rho;
    u2(right) = rhoV;
    u3(right) = rhoU;
    u4(right) = rhoE;
end

end
