function [F1, F2, F3] = SolveFluxRoe(U1L, U2L, U3L, U1R, U2R, U3R, Gamma)

% Переход к консервативным переменным слева и справа:
rhoL = U1L;
uL = U2L./rhoL;
pL = (Gamma - 1) * (U3L - 0.5 * (uL.^2).*rhoL);
cL = sqrt((Gamma * pL./rhoL));
HL = cL.^2 / (Gamma - 1) + uL.^2 / 2;
rhoR = U1R;
uR = U2R./rhoR;
pR = (Gamma - 1) * (U3R - 0.5 * (uR.^2).*rhoR);
cR = sqrt((Gamma * pR./rhoR));
HR = cR.^2 / (Gamma - 1) + uR.^2 / 2;

% Осредненные параметры (Роэ):
u_bar = (sqrt(rhoL).*uL + sqrt(rhoR).*uR)./(sqrt(rhoL) + sqrt(rhoR));
H_bar = (sqrt(rhoL).*HL + sqrt(rhoR).*HR)./(sqrt(rhoL) + sqrt(rhoR));
c_bar = sqrt((Gamma - 1) * (H_bar - 0.5 * u_bar.^2));

% dU = U(R) - U(L) на границах ячеек {x_{i}, x_{i+1}, ...}:
dU1 = U1R - U1L;
dU2 = U2R - U2L;
dU3 = U3R - U3L;

% Потоковые члены слева и справа F_{i+1/2}:
F1L = rhoL.*uL; F2L = rhoL.*uL.^2 + pL; F3L = rhoL.*uL.*HL;
FL = [F1L; F2L; F3L];
F1R = rhoR.*uR; F2R = rhoR.*uR.^2 + pR; F3R = rhoR.*uR.*HR;
FR = [F1R; F2R; F3R];

% Правые собственные вектора:
r1 = [1*ones(1,length(u_bar)); u_bar-c_bar; H_bar-u_bar.*c_bar];
r2 = [1*ones(1,length(u_bar)); u_bar; 1/2*u_bar.^2];
r3 = [1*ones(1,length(u_bar)); u_bar+c_bar; H_bar+u_bar.*c_bar];

% Характеристические приращения:
dw2 = (Gamma - 1)./(c_bar.^2).*(dU1.*(H_bar - u_bar.^2) + u_bar.*dU2 - dU3);
dw1 = 1./(2 * c_bar).*(dU1.*(u_bar + c_bar) - dU2 - c_bar.*dw2);
dw3 = dU1 - (dw1 + dw2);

% Собственные числа:
lambda1 = abs(u_bar - c_bar);
lambda2 = abs(u_bar);
lambda3 = abs(u_bar + c_bar);
lambda1L = uL - cL;
lambda3L = uL + cL;
lambda1R = uR - cR;
lambda3R = uR + cR;

% Искусственная вязкость вблизи звуковой точки:
del = 0.0001;
boolv1 = (lambda1R - lambda1L) > 0;
eps = 2.0 * (lambda1R - lambda1L).*boolv1;
boolv1 = lambda1 < eps;
lambda1 = 0.5 * (eps + (lambda1.^2) / (eps + del)).*boolv1 + lambda1.*(1 - boolv1);
boolv2 = (lambda3R - lambda3L) > 0;
eps = 2.0 * (lambda3R - lambda3L).*boolv2;
boolv2 = lambda3 < eps;
lambda3 = 0.5 * (eps + (lambda3.^2) / (eps + del)).*boolv2 + lambda3.*(1 - boolv2);

% Приведение данных:
dw1 = repmat(dw1, 3, 1);
dw2 = repmat(dw2, 3, 1);
dw3 = repmat(dw3, 3, 1);
lambda1 = repmat(lambda1, 3, 1);
lambda2 = repmat(lambda2, 3, 1);
lambda3 = repmat(lambda3, 3, 1);

% Вычисление потоков (Роэ):
Flux = 0.5 * (FL + FR) - 0.5 * (dw1.*lambda1.*r1 + dw2.*lambda2.*r2 + dw3.*lambda3.*r3);
F1 = Flux(1,:);
F2 = Flux(2,:);
F3 = Flux(3,:);

end
