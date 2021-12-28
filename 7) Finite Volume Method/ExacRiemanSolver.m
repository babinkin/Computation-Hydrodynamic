% Функция расчета параметров задачи Римана для уравнений Эйлера (задача о
% распаде произвольного разрыва) точным подходом
%
%
% используя аналит соотношения можно составить неленейное уравнение для
% определения давления на контактном разрыве, решив его итерационным методом
% Ньютона (метод касательных) можно пересчитать остальные параметры во всех
% областях.
% 
% Входные параметры:
%
% left_density              - плотность слева от разрыва
% left_velocity             - скорость --//--
% left_pressure             - давление --//--                                            
% right_density             - плотность справа от разрыва
% right_velocity            - скорость --//--
% right_pressure            - давление --//--
% Gamma                     - показатель адиабаты
% lambda                    - автомодельная переменная lambda = x/t
% MaxIteration              - максимальное количество итераций
% TOL                       - константа условия сходимости
%
%
% Выходные параметры:
%
% ret_density               - плотность на линии x/t = lambda
% ret_velocity              - скорость --//--
% ret_pressure              - давление --//--
% is_left_of_contact        - признак слева ли от контактного разрыва
% iteration                 - выполненое количество итераций
%
% =========================================================================
function [ret_density,ret_velocity,ret_pressure,is_left_of_contact,iteration] = ExacRiemanSolver(left_density,left_velocity,left_pressure,...
    right_density,right_velocity,right_pressure,Gamma, lambda,MaxIteration,TOL)

iteration = 0;
% Все формулы взяты из Торо
%==========================================================================
% solve Sound Velosity for Left and Right
left_soundspeed=sqrt ( Gamma*left_pressure/left_density );
right_soundspeed=sqrt( Gamma*right_pressure/right_density );
%  проверяем образуется ли вакуум (образуется редко)
left_vacuum_front_speed = left_velocity + 2.0 * left_soundspeed / ( Gamma - 1.0 );
right_vacuum_front_speed = right_velocity - 2.0 * right_soundspeed / ( Gamma - 1.0 );
critical_speed =  left_vacuum_front_speed - right_vacuum_front_speed;
%---------------------------------------------------------
if ( critical_speed < 0.0 ) % образуется зона вукуума
    left_head_speed = left_velocity - left_soundspeed;
	left_tail_speed = left_vacuum_front_speed;
	right_head_speed = right_velocity + right_soundspeed;
	right_tail_speed = right_vacuum_front_speed;
	%-----------------------
    is_left_of_contact = lambda < left_tail_speed;
	if ( is_left_of_contact ) % определяем где находится искомая линия lambda слева ли от контактного разрыва
        if ( lambda < left_head_speed ) 
            ret_density  = left_density;
			ret_velocity = left_velocity;
			ret_pressure = left_pressure;
		else 
			% left_rarefaction (4.56)
			temp1 = 2.0/ ( Gamma+1.0 ) + ( Gamma-1.0 ) / ( Gamma+1.0 ) /left_soundspeed *(left_velocity - lambda);
			ret_density = left_density * temp1^( 2.0/ ( Gamma-1.0 ) );
			ret_pressure = left_pressure *  temp1^( 2.0*Gamma/ ( Gamma-1.0 ) );
			ret_velocity = 2.0/ ( Gamma+1.0 ) * ( left_soundspeed + ( Gamma-1.0 ) /2.0*left_velocity + lambda );
        end
    else
		if ( lambda > right_tail_speed ) 
			if ( lambda > right_head_speed )
				ret_density  = right_density;
				ret_velocity = right_velocity;
				ret_pressure = right_pressure;
            else
				%right_rarefaction (4.63)
				temp1 =2.0/ ( Gamma+1.0 ) - ( Gamma-1.0 ) / ( Gamma+1.0 ) /right_soundspeed *(right_velocity - lambda);
				ret_density = right_density * temp1^( 2.0/ ( Gamma-1.0 ) );
				ret_pressure = right_pressure *  temp1^(2.0*Gamma/ ( Gamma-1.0 ) );
				ret_velocity = 2.0/ ( Gamma+1.0 ) * ( -right_soundspeed + ( Gamma-1.0 ) /2.0*right_velocity + lambda);
            end
        else
			% u resides inside vaccum
			ret_density=0.0;
            ret_velocity=0.0;
            ret_pressure = 0.0;
        end
    end
%---------------------------------------------------------
else
	%Решение итерационным способом
	%Начальное приближение
	%-------------------------------
	% PVRS (9.20)
	p_star=0.5* ( left_pressure+right_pressure ) +0.125* ( left_velocity-right_velocity ) * ( left_density+right_density ) * ( left_soundspeed+right_soundspeed );
	p_star=max(p_star,TOL);
	%----------
	pMin=min(left_pressure,right_pressure);
	pMax=max(left_pressure,right_pressure);
	if ( p_star>pMax ) 
		% TSRS (9.42)
		temp1=sqrt ( ( 2.0/ ( Gamma+1.0 ) /left_density ) / ( p_star+ ( Gamma-1.0 ) / ( Gamma+1.0 ) *left_pressure ) );
		temp2=sqrt ( ( 2.0/ ( Gamma+1.0 ) /right_density ) / ( p_star+ ( Gamma-1.0 ) / ( Gamma+1.0 ) *right_pressure ) );
		p_star= ( temp1*left_pressure+temp2*right_pressure+ ( left_velocity-right_velocity ) ) / ( temp1+temp2 );
		p_star=max(p_star,TOL);
    else
        if ( p_star<pMin )
			%  TRRS (9.32)
			temp1= ( Gamma-1.0 ) / ( 2.0*Gamma );
			p_star= ( ( left_soundspeed+right_soundspeed+0.5* ( Gamma-1.0 ) * ( left_velocity-right_velocity ) ) / ( left_soundspeed/left_pressure^temp1 +...
                right_soundspeed/right_pressure^temp1 ) )^( 1.0/temp1 );
        else
            % p_star = PVRS
        end
    end
	%---------------------------------------------------------------
	%main loop for Pressure
	for iteration=1:MaxIteration
        %LEFT
		% solve temp value
		temp1=sqrt ( ( 2.0/ ( Gamma+1.0 ) /left_density ) / ( p_star+ ( Gamma-1.0 ) / ( Gamma+1.0 ) *left_pressure ) );
		% ---Solve function 4.6 4.7
        if p_star<=left_pressure
            f1=2.0/ ( Gamma-1.0 ) *left_soundspeed* (    ( p_star/left_pressure)^( ( Gamma-1.0 ) / ( 2.0*Gamma ) )   - 1.0 ) ;
        else
            f1= ( p_star-left_pressure ) *temp1;
        end
		%Solve derivates 4.37
        if p_star<=left_pressure
            f_d= ( p_star/left_pressure)^(- ( Gamma+1.0 ) / ( 2.0*Gamma ) ) / ( left_density*left_soundspeed );
        else
            f_d=temp1* ( 1.0-0.5* ( p_star-left_pressure ) / ( p_star+ ( Gamma-1.0 ) / ( Gamma+1.0 ) *left_pressure ) );
        end
		%RIGHT
		temp1=sqrt ( ( 2.0/ ( Gamma+1.0 ) /right_density ) / ( p_star+ ( Gamma-1.0 ) / ( Gamma+1.0 ) *right_pressure ) );
		%---Solve function 4.6 4.7
        if p_star<=right_pressure
            f2=  2.0/ ( Gamma-1.0 ) *right_soundspeed* (    ( p_star/right_pressure)^( ( Gamma-1.0 ) / ( 2.0*Gamma ) )   -1.0    ) ;
        else
            f2= ( p_star-right_pressure ) *temp1;
        end
		%Solve derivates
        if p_star<=right_pressure
            f_d=f_d+  ( p_star/right_pressure)^(- ( Gamma+1.0 ) / ( 2.0*Gamma ) ) / ( right_density*right_soundspeed ) ;
        else
            f_d=f_d+  temp1* ( 1.0-0.5* ( p_star- right_pressure ) / ( p_star+ ( Gamma-1.0 ) / ( Gamma+1.0 ) *right_pressure ) );
        end
		%------------------------------------------------------------------
		p_New=p_star- ( f1+f2- ( left_velocity-right_velocity ) ) /f_d;
        if ( abs ( p_New - p_star ) / ( 0.5 * abs ( p_New + p_star ) ) < TOL ) 
            break
        end
		p_star=p_New;
    end
	%-----------------------------------------------------------------------------------------------------
    %-----------------------------------------------------------------------------------------------------
	% calculate star speed */
	star_speed=0.5* ( left_velocity + right_velocity ) +0.5* ( f2-f1 );
	% calculate other star values */
	% Left
	if ( p_star>=left_pressure ) 
        % SHOCK
		left_star_density = left_density * ( p_star / left_pressure + ( Gamma-1.0 ) / ( Gamma+1.0 ) ) /...
                ( ( Gamma-1.0 ) / ( Gamma+1.0 ) * p_star / left_pressure + 1.0 );
		left_tail_speed = left_velocity -left_soundspeed * sqrt ( ( Gamma+1.0 ) / ( 2.0*Gamma ) * p_star/left_pressure + ...
                ( Gamma-1.0 ) / ( 2.0*Gamma ) );
        left_head_speed = left_tail_speed;
	else  % left_wave_ == kRarefaction
		left_star_density = left_density * ( p_star / left_pressure)^(1.0/Gamma );
		left_head_speed = left_velocity - left_soundspeed;
		left_tail_speed = star_speed - sqrt ( Gamma*p_star/left_star_density );
    end
	% Right
	if ( p_star>=right_pressure ) 
		right_star_density = right_density *...
			                     ( p_star / right_pressure + ( Gamma-1.0 ) / ( Gamma+1.0 ) ) /...
			                     ( ( Gamma-1.0 ) / ( Gamma+1.0 ) * p_star / right_pressure + 1.0 );
		right_tail_speed = right_velocity +...
                right_soundspeed * sqrt ( ( Gamma+1.0 ) / ( 2.0*Gamma ) * p_star/right_pressure +...
                ( Gamma-1.0 ) / ( 2.0*Gamma ) );
        right_head_speed = right_tail_speed;
    else  % right_wave_ == kRarefaction
		right_star_density = right_density *  ( p_star / right_pressure)^ ( 1.0/Gamma );
		right_head_speed = right_velocity + right_soundspeed;
		right_tail_speed = star_speed + sqrt ( Gamma*p_star/right_star_density );
    end
    %-----------------------------------------------------------------------------------------------------
    % LEFT
    is_left_of_contact = lambda  < star_speed;
	if ( is_left_of_contact ) % the u is left of contact discontinuity
        if ( p_star>=left_pressure )  % the left wave is a shock
            if ( lambda < left_head_speed )  % the u is before the shock
                ret_density  = left_density;
				ret_velocity = left_velocity;
				ret_pressure = left_pressure;
			else  % the u is behind the shock
				ret_density  = left_star_density;
				ret_velocity = star_speed;
				ret_pressure = p_star;
            end
		else  % the left wave is a rarefaction
			if ( lambda < left_head_speed )  % the u is before the rarefaction
				ret_density  = left_density;
				ret_velocity = left_velocity;
				ret_pressure = left_pressure;
            else
				if ( lambda < left_tail_speed )  % the u is inside the rarefaction
                    % left_rarefaction (4.56)
					temp1 = 2.0/ ( Gamma+1.0 ) + ( Gamma-1.0 ) / ( Gamma+1.0 ) /left_soundspeed *(left_velocity - lambda);
					ret_density = left_density *  temp1^( 2.0/ ( Gamma-1.0 ) );
					ret_pressure = left_pressure * temp1^( 2.0*Gamma/ ( Gamma-1.0 ) );
					ret_velocity = 2.0/ ( Gamma+1.0 ) * ( left_soundspeed + ( Gamma-1.0 ) /2.0*left_velocity + lambda);
				else  % the u is after the rarefaction
					ret_density  = left_star_density;
					ret_velocity = star_speed;
					ret_pressure = p_star;
                end
            end
        end
    else % the queried u is right of contact discontinuity
        %------------------------------------------------------------------------
		if ( p_star>=right_pressure )  % the right wave is a shock
			if ( lambda > right_head_speed )  % the u is before the shock
				ret_density  = right_density;
				ret_velocity = right_velocity;
				ret_pressure = right_pressure;
			else  % the u is behind the shock
				ret_density  = right_star_density;
				ret_velocity = star_speed;
				ret_pressure = p_star;
            end
        else  % the right wave is a rarefaction
			if ( lambda > right_head_speed )  % the u is before the rarefaction
				ret_density  = right_density;
				ret_velocity = right_velocity;
				ret_pressure = right_pressure;
			else 
				if ( lambda > right_tail_speed )  % the u is inside the rarefaction
					%right_rarefaction (4.63)
					temp1 =2.0/ ( Gamma+1.0 ) - ( Gamma-1.0 ) / ( Gamma+1.0 ) /right_soundspeed *(right_velocity - lambda);
					ret_density = right_density *  temp1^( 2.0/ ( Gamma-1.0 ) );
					ret_pressure = right_pressure * temp1^( 2.0*Gamma/ ( Gamma-1.0 ) );
					ret_velocity = 2.0/ ( Gamma+1.0 ) * ( -right_soundspeed + ( Gamma-1.0 ) /2.0*right_velocity + lambda);
				else  % the u is after the rarefaction
					ret_density  = right_star_density;
					ret_velocity = star_speed;
					ret_pressure = p_star;
                end
            end
        end        
    end
end

end