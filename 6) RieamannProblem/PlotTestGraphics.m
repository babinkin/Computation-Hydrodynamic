function PlotTestGraphics(~, U1new, U2new, U3new, rhoL, uL, pL, rhoR, uR, pR, x_i, t, x_dis, Gamma, R, TextSize, test, PlotType)

rho_plot = U1new;
u_plot = U2new./rho_plot;
p_plot = (Gamma - 1) * (U3new - 0.5 * (u_plot.^2).*rho_plot);
T_plot = p_plot./(R * rho_plot);
int_energ_plot = (Gamma) / (Gamma - 1) * p_plot./rho_plot;
rho_exac = 0;
u_exac = 0;
p_exac = 0;
int_energ_exac = 0;

if (test <= 0 || (test == 6))

else
    
    [rho_exac, u_exac, p_exac] = ExacRiemanProblemSolution(rhoL, uL, pL, rhoR, uR, ...
        pR, x_i, t(end), x_dis, Gamma);
    int_energ_exac = (Gamma) / (Gamma - 1) * p_exac./rho_exac;
    T_exac = p_exac./(R * rho_exac);

end

switch (PlotType)
    case {'Rho'}
       
        var_plot = rho_plot;
        exac_plot = rho_exac;
    
    case {'Vel'}
       
        var_plot = u_plot;
        exac_plot = u_exac;
   
    case {'Press'}
        
        var_plot = p_plot;
        exac_plot = p_exac;
    
    case {'IntEnerg'}
       
        var_plot = int_energ_plot;
        exac_plot = int_energ_exac;
    
    case {'Temp'}
        
        var_plot = T_plot;
        exac_plot = T_exac;
    
    case {'RhoV'}
        
        var_plot = u_plot.*rho_plot;
        exac_plot = u_exac.*rho_exac;

    otherwise
        
        var_plot = rho_plot;
        exac_plot = rho_exac;

end

PlotVar(var_plot, exac_plot, x_i, t, TextSize, test);
pause(0.05);

end

function PlotVar(var_plot, exac_plot, x_i, t, TextSize, test)

figure(1);

if ((test <= 0) || (test == 6))
    
    plot(x_i, var_plot(2:end-1), 'k.');

    if (test == 6)
       
        xlim([0.5, 0.9]);

    end
    
else
    plot(x_i, var_plot(2:end-1), '-r.', x_i, exac_plot,'-k');
    minimum = min(var_plot);
    maximum = max(var_plot);
    ylim([minimum * (1 - 0.1 * sign(minimum)), maximum * (1 + 0.1 * sign(maximum))]);
    legend('Solver', 'Analytical', 'Location', 'southwest', 'FontSize', TextSize-10);

end

grid on;
title(['Time = ', num2str(t(end)), ' '], 'FontSize', TextSize);
xlabel('x', 'FontSize', TextSize);
ylabel('Value', 'FontSize', TextSize);
set(gca, 'FontSize', 18);

end







