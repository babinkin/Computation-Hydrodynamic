% Функция строит точное решение для линейного уравнения переноса со
% ступенькой с учетом переодичных ГУ FIXME
%
% Входные параметры:
% x         - массив координат узлов
% t         - время
% Uinit 	- высота ступеньки
% xStartStep - координата от которой начинается ступенька в нач момент врем
% Lstep     - ширина ступеньки
%
%
%
% Выходные параметры:
%
% U         - точное решение на заданной сетке узлов x
%
% =========================================================================
function U = U_Exac_Step(x,t,V,Uinit,xStartStep,Lstep,xLeft,xRight)
% len = length(x);
L = xRight - xLeft;
% xRight % (V*t)
n  = floor((xStartStep + V*t)/xRight);
x1 = xStartStep + V*t - n*L 
x2 = xStartStep + Lstep + V*t - n*L
% % ---



if ( (x2 > xRight)  )
    x2 = xLeft + (x2 -xRight )
    U = f2(x,x1,x2,Uinit,xLeft,xRight);
else
    U = f1(x,x1,x2,Uinit);
end


end


function f1 = f1(x,x1,x2,Uinit)
len = length(x);
x_c = 0.5*(x(1:end-1) + x(2:end));
ind =  find(x_c>=x1 & x_c<=x2)
f1 =  zeros(1,len);
f1(ind) = Uinit;
f1(ind+1) = Uinit;
end

function f2 = f2(x,x1,x2,Uinit,xLeft,xRight)
len = length(x);
f2 =  zeros(1,len);
x_c = 0.5*(x(1:end-1) + x(2:end));
ind =  find( (x_c>=x1 & x_c<=xRight) );
f2(ind) = Uinit;
f2(ind+1) = Uinit;
% -- 
ind =  find(  (x_c>=xLeft & x_c<=x2) );
f2(ind) = Uinit;
end



















