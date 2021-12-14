function [rho, u, p] = ExacRiemanProblemSolution(rhoL, uL, pL, rhoR, uR, pR, x_i, t, x0, Gamma)

MaxIteration = 20;
TOL = 1e-9;
len = length(x_i);
rho = zeros(1, len);
u = zeros(1, len);
p = zeros(1, len);

for i = 1:len
   
    lambda = (x_i(i) - x0) / t;
    [rho(i), u(i), p(i), ~, ~] = ExacRiemanSolver(rhoL, uL, pL,rhoR, uR, pR, Gamma, lambda, MaxIteration, TOL);

end

end




















