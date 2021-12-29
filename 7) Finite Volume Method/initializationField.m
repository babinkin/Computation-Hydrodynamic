% функция производит инициализацию полей
%
%
%
%
% =========================================================================
function [u1,u2,u3,u4] = initializationField(pres,vx,vy,T,k,Rgas,numCell)

u1 = ones(numCell,1)*pres/(Rgas*T);
u2 = u1*vx;
u3 = u2*vy;
u4 = pres/(k-1) + 0.5*u1*(vx*vx + vy*vy);





end