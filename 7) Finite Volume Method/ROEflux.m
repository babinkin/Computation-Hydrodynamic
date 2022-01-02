function Roe = ROEflux(left_density,left_velocity,left_E,left_pressure,right_density,right_velocity,right_E,right_pressure,gamma)
    % Roe flux function
    %
    
    % Left state
    rL = left_density;
    uL = left_velocity;
    EL = left_E./rL;
    pL = left_pressure;
    aL = sqrt(gamma*pL/rL);
    HL = ( left_E + pL )./rL;
    
    % Right state
    rR = right_density;
    uR = right_velocity;
    ER = right_E./rR;
    pR = right_pressure;
    aR = sqrt(gamma*pR/rR);
    HR = ( right_E + pR )./rR;
    
    % First compute the Roe Averages
    RT = sqrt(rR/rL);
    r = RT*rL;
    u = (uL+RT*uR)/(1+RT);
    H = (HL+RT*HR)/(1+RT);
    a = sqrt( (gamma-1)*(H-u*u/2) );
    
    % Differences in primitive variables.
    dr = rR - rL;
    du = uR - uL;
    dP = pR - pL;
    
    % Wave strength (Characteristic Variables).
    dV = [(dP-r*a*du)/(2*a^2); -( dP/(a^2)-dr); (dP+r*a*du)/(2*a^2)];
    
    % Absolute values of the wave speeds (Eigenvalues)
    ws = [ abs(u-a); abs( u ); abs(u+a) ];

    
    % Right eigenvectors
    R = [  1  ,  1  ,  1  ;
         u-a ,  u  , u+a ;
        H-u*a,u^2/2,H+u*a];
   
    % Compute the average flux.
    FL=[rL.*uL; rL.*uL.^2+pL; uL.*(rL.*EL+pL)];
    FR=[rR.*uR; rR.*uR.^2+pR; uR.*(rR.*ER+pR)];

    % Add the matrix dissipation term to complete the Roe flux.
    Roe = ( FL + FR  - R*(ws.*dV))/2;
end
