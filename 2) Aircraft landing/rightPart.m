% Функция правых частей уравнения движения ЛА в поле земли
function F = rightPart(t,Y,rho0,g0,rz,sigma,K)

% достаем переменные
V = Y(1);
theta = Y(2);
x = Y(3);
h = Y(4);

% Вычисляем силы 
rho = rho0*exp(-h/7800);
g = g0*(rz./(rz+h)).^2;
q = 0.5*rho.*(V.^2);

% вычисляем правые части
F    = [ -sigma * rho * V^2 / 2 - g * sin(theta);
     (1 / V) * (sigma * K * rho * V^2 / 2 - g * cos(theta) + ...
        V^2 / (rz + h) * cos(theta));
     V * cos(theta);
     V * sin(theta);
     ];

end
