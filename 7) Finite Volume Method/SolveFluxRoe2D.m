function [F1,F2,F3,F4] = SolveFluxRoe2D(U1L,U2L,U3L,U4L,U1R,U2R,U3R,U4R,Gamma)

%
% U1L,U2L,U3L           - вектор консервативных переменных слева
% U1R,U2R,U3R           - вектор консервативных переменных справа
% 
%
%
% Выходные параметры
% F1,F2,F3              - вектор потока на грани КО
% iteration             - выполненое количество итераций
%
% =========================================================================


% ---
% переводим в физические переменные
left_density=U1L;
right_density=U1R;
left_velocity=U2L./left_density;
right_velocity=U2R./right_density;
v3L = U3L./left_density;
v3R = U3R./right_density;
left_pressure=(Gamma-1)*(U4L-0.5*(left_velocity.^2+v3L.^2).*left_density);
right_pressure=(Gamma-1)*(U4R-0.5*(right_velocity.^2 + v3R.^2).*right_density);
% ----------------
MaxIteration=20; % макс число итераций
TOL=1e-8;
lambda = 0; % линия на грани КО
% ---
len = length(U1L);

for i=1:len
Roe(:,i) = ROEflux(left_density(i),left_velocity(i),U4L(i),left_pressure(i),right_density(i),right_velocity(i),U4R(i),right_pressure(i),Gamma);
end


% расчитываем потоки
F1=Roe(1,:);
F2=Roe(2,:);
ind = (F1>=0);
v3 = (v3L').*ind + (v3R').*(~ind);
F3 = F1.*v3;
F4=Roe(3,:);

end
