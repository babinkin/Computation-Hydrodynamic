% максимально возможная скорость истечения
function f=Vmax_ot_T0(T0,R,k)
f=(2*k/(k-1)*R*T0).^(0.5);