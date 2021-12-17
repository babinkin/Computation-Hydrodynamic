function [U1new, U2new, U3new] = SolveEvolutionExplFirstOrder(F1, F2, F3, U1old, U2old, U3old, nu, dt, rf_dr, rc_int, gamma)
len = length(F1);
len_dr= length(rf_dr);
% Координаты граней:

if len_dr < 2
   
    rf = ((1:len) - 1) * rf_dr;
    rc_int = rf(1:end-1) + rf_dr / 2;

else
    
    rf = rf_dr;

end

U1new = zeros(1,len+1);
U2new = zeros(1,len+1);
U3new = zeros(1,len+1);

temp = dt * nu./(rf(2:end).^nu - rf(1:end-1).^nu);

U1new(2:end-1) = U1old(2:end-1) - temp.*(F1(2:end).*(rf(2:end).^(nu-1)) - F1(1:end-1).*(rf(1:end-1).^(nu-1)));
U2new(2:end-1) = U2old(2:end-1) - temp.*(F2(2:end).*(rf(2:end).^(nu-1)) - F2(1:end-1).*(rf(1:end-1).^(nu-1))) + (nu-1)*(gamma-1)*temp.*(U3old(2:end-1) - rf(1:end-1) ).*(rc_int.^(nu-2));
U3new(2:end-1) = U3old(2:end-1) - temp.*(F3(2:end).*(rf(2:end).^(nu-1)) - F3(1:end-1).*(rf(1:end-1).^(nu-1)));
end





