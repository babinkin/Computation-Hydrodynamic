% функция решает нелинейное уравнеие для расходной функции
% 
% M_start - начальное приближение
%

function M = SolveQFun(q, k, M0) 

    M = M0 - f(q, k, M0) / df(k, M0);
    
    while abs(M - M0) > 1e-6
        
        M0 = M;
        M = M0 - f(q, k, M0) / df(k, M0);
        
    end

end
