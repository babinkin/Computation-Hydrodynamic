% Скорость газа от числа Маха
function V=V_ot_M_T0_R_k(M,T0,R,k)
V=Vmax_ot_T0(T0,R,k)*V_na_Vmax_ot_M(M,k);
end