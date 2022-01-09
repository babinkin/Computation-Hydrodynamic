% Функция расчета шага по времени для явной схемы эйлера
%
%
%
% -------------------------------------------------------------------------
function dt = SolveTimeStepExplicitEiler(U1,U2,U3,dx,Courant,gamma)
rho = U1;
u = U2./rho;
p = (gamma-1)*(U3 - 0.5*(u.^2).*rho);
c = sqrt(gamma*p./rho);
dt = min( Courant*dx./(abs(u)+c) );
end