%V/Vmax
function f=V_na_Vmax_ot_M(M,k)
f=(((k-1)/2.*M.^2)./(1+(k-1)/2.*M.^2)).^(0.5);