classdef trebfunctions 
methods(Static) 
function [T_p] = T_p(Yct, Ycdt, thSt, thSdt, l_a, l_s) 
	T_p = (2.25*(l_a*(-Ycdt*Yct*(-Yct^2/4 - 1.5)/2 + Ycdt*Yct/4)*sin(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) + l_s*thSdt*cos(thSt))^2 + 2.25*(-l_a*(-Ycdt*Yct*(-Yct^2/4 - 1.5)/2 + Ycdt*Yct/4)*cos(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - l_s*thSdt*sin(thSt))^2);
end


function [Ycam] = Ycam(Yct, Ycdt, thSt, thSdt, l_a, l_s) 
	Ycam = (0);
end


function [thSdd] = thSdd(Yct, Ycdt, thSt, thSdt, l_a, l_s) 
	thSdd = ((4.5*l_s*(l_a*(-Ycdt*Yct*(-Yct^2/4 - 1.5)/2 + Ycdt*Yct/4)^2*sin(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - l_a*(Ycdt^2*Yct^2/4 - Ycdt^2*(-Yct^2/4 - 1.5)/2 + Ycdt^2/4)*cos(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - l_a*(-Ycdt*Yct*(-Yct^2/4 - 1.5)/2 + Ycdt*Yct/4)*(Ycdt*Yct*(-Yct^2/4 - 1.5)/2 - Ycdt*Yct/4)*cos(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/(Yct^2/4 + (-Yct^2/4 - 1.5)^2)^(3/2) - l_s*thSdt^2*cos(thSt))*sin(thSt) - 4.5*l_s*(l_a*(-Ycdt*Yct*(-Yct^2/4 - 1.5)/2 + Ycdt*Yct/4)^2*cos(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/(Yct^2/4 + (-Yct^2/4 - 1.5)^2) + l_a*(Ycdt^2*Yct^2/4 - Ycdt^2*(-Yct^2/4 - 1.5)/2 + Ycdt^2/4)*sin(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) + l_a*(-Ycdt*Yct*(-Yct^2/4 - 1.5)/2 + Ycdt*Yct/4)*(Ycdt*Yct*(-Yct^2/4 - 1.5)/2 - Ycdt*Yct/4)*sin(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/(Yct^2/4 + (-Yct^2/4 - 1.5)^2)^(3/2) - l_s*thSdt^2*sin(thSt))*cos(thSt) - 44.1*l_s*cos(thSt) - 2.25*(2*l_a*l_s*(-Yct*(-Yct^2/4 - 1.5)/2 + Yct/4)*sin(thSt)*cos(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) + 2*l_a*l_s*(-Yct*(-Yct^2/4 - 1.5)/2 + Yct/4)*sin(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)*cos(thSt)/sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2))*(-l_a^2*(-Yct*(-Yct^2/4 - 1.5) + Yct/2)*(Ycdt^2*Yct^2/4 - Ycdt^2*(-Yct^2/4 - 1.5)/2 + Ycdt^2/4)/(2*(Yct^2/4 + (-Yct^2/4 - 1.5)^2)) - l_a^2*(-Yct*(-Yct^2/4 - 1.5) + Yct/2)*(-Ycdt*Yct*(-Yct^2/4 - 1.5)/2 + Ycdt*Yct/4)*(Ycdt*Yct*(-Yct^2/4 - 1.5) - Ycdt*Yct/2)/(2*(Yct^2/4 + (-Yct^2/4 - 1.5)^2)^2) + l_a^2*(Yct*(-Yct^2/4 - 1.5) - Yct/2)*(-Ycdt*Yct*(-Yct^2/4 - 1.5)/2 + Ycdt*Yct/4)^2/(2*(Yct^2/4 + (-Yct^2/4 - 1.5)^2)^2) - 58.8*l_a^2*(-Yct*(-Yct^2/4 - 1.5)/2 + Yct/4)*sin(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 4.5*l_a*(-Yct*(-Yct^2/4 - 1.5)/2 + Yct/4)*(-Ycdt*Yct*(-Yct^2/4 - 1.5)/2 + Ycdt*Yct/4)*(l_a*(-Ycdt*Yct*(-Yct^2/4 - 1.5)/2 + Ycdt*Yct/4)*sin(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) + l_s*thSdt*cos(thSt))*cos(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 4.5*l_a*(-Yct*(-Yct^2/4 - 1.5)/2 + Yct/4)*(-Ycdt*Yct*(-Yct^2/4 - 1.5)/2 + Ycdt*Yct/4)*(-l_a*(-Ycdt*Yct*(-Yct^2/4 - 1.5)/2 + Ycdt*Yct/4)*cos(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - l_s*thSdt*sin(thSt))*sin(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/(Yct^2/4 + (-Yct^2/4 - 1.5)^2) + 4.5*l_a*(-Yct*(-Yct^2/4 - 1.5)/2 + Yct/4)*(l_a*(-Ycdt*Yct*(-Yct^2/4 - 1.5)/2 + Ycdt*Yct/4)^2*sin(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - l_a*(Ycdt^2*Yct^2/4 - Ycdt^2*(-Yct^2/4 - 1.5)/2 + Ycdt^2/4)*cos(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - l_a*(-Ycdt*Yct*(-Yct^2/4 - 1.5)/2 + Ycdt*Yct/4)*(Ycdt*Yct*(-Yct^2/4 - 1.5)/2 - Ycdt*Yct/4)*cos(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/(Yct^2/4 + (-Yct^2/4 - 1.5)^2)^(3/2) - l_s*thSdt^2*cos(thSt))*cos(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 4.5*l_a*(-Yct*(-Yct^2/4 - 1.5)/2 + Yct/4)*(l_a*(-Ycdt*Yct*(-Yct^2/4 - 1.5)/2 + Ycdt*Yct/4)^2*cos(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/(Yct^2/4 + (-Yct^2/4 - 1.5)^2) + l_a*(Ycdt^2*Yct^2/4 - Ycdt^2*(-Yct^2/4 - 1.5)/2 + Ycdt^2/4)*sin(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) + l_a*(-Ycdt*Yct*(-Yct^2/4 - 1.5)/2 + Ycdt*Yct/4)*(Ycdt*Yct*(-Yct^2/4 - 1.5)/2 - Ycdt*Yct/4)*sin(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/(Yct^2/4 + (-Yct^2/4 - 1.5)^2)^(3/2) - l_s*thSdt^2*sin(thSt))*sin(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 44.1*l_a*(-Yct*(-Yct^2/4 - 1.5)/2 + Yct/4)*sin(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 4.5*l_a*(l_a*(-Ycdt*Yct*(-Yct^2/4 - 1.5)/2 + Ycdt*Yct/4)*sin(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) + l_s*thSdt*cos(thSt))*(Ycdt*Yct^2/4 - Ycdt*(-Yct^2/4 - 1.5)/2 + Ycdt/4)*sin(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) + 4.5*l_a*(-l_a*(-Ycdt*Yct*(-Yct^2/4 - 1.5)/2 + Ycdt*Yct/4)*cos(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - l_s*thSdt*sin(thSt))*(Ycdt*Yct^2/4 - Ycdt*(-Yct^2/4 - 1.5)/2 + Ycdt/4)*cos(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 4.5*l_a*(-Yct*(-Yct^2/4 - 1.5)/2 + Yct/4)*(Ycdt*Yct*(-Yct^2/4 - 1.5)/2 - Ycdt*Yct/4)*(l_a*(-Ycdt*Yct*(-Yct^2/4 - 1.5)/2 + Ycdt*Yct/4)*sin(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) + l_s*thSdt*cos(thSt))*sin(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/(Yct^2/4 + (-Yct^2/4 - 1.5)^2)^(3/2) + 4.5*l_a*(-Yct*(-Yct^2/4 - 1.5)/2 + Yct/4)*(Ycdt*Yct*(-Yct^2/4 - 1.5)/2 - Ycdt*Yct/4)*(-l_a*(-Ycdt*Yct*(-Yct^2/4 - 1.5)/2 + Ycdt*Yct/4)*cos(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - l_s*thSdt*sin(thSt))*cos(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/(Yct^2/4 + (-Yct^2/4 - 1.5)^2)^(3/2) + 2.25*(l_a*(-Ycdt*Yct*(-Yct^2/4 - 1.5)/2 + Ycdt*Yct/4)*sin(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) + l_s*thSdt*cos(thSt))*(2*l_a*(-Yct*(-Yct^2/4 - 1.5)/2 + Yct/4)*(-Ycdt*Yct*(-Yct^2/4 - 1.5)/2 + Ycdt*Yct/4)*cos(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/(Yct^2/4 + (-Yct^2/4 - 1.5)^2) + 2*l_a*(Ycdt*Yct^2/4 - Ycdt*(-Yct^2/4 - 1.5)/2 + Ycdt/4)*sin(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) + 2*l_a*(Yct*(-Yct^2/4 - 1.5)/2 - Yct/4)*(-Ycdt*Yct*(-Yct^2/4 - 1.5)/2 + Ycdt*Yct/4)*sin(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/(Yct^2/4 + (-Yct^2/4 - 1.5)^2)^(3/2)) + 2.25*(-l_a*(-Ycdt*Yct*(-Yct^2/4 - 1.5)/2 + Ycdt*Yct/4)*cos(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - l_s*thSdt*sin(thSt))*(2*l_a*(-Yct*(-Yct^2/4 - 1.5)/2 + Yct/4)*(-Ycdt*Yct*(-Yct^2/4 - 1.5)/2 + Ycdt*Yct/4)*sin(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 2*l_a*(Ycdt*Yct^2/4 - Ycdt*(-Yct^2/4 - 1.5)/2 + Ycdt/4)*cos(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 2*l_a*(Yct*(-Yct^2/4 - 1.5)/2 - Yct/4)*(-Ycdt*Yct*(-Yct^2/4 - 1.5)/2 + Ycdt*Yct/4)*cos(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/(Yct^2/4 + (-Yct^2/4 - 1.5)^2)^(3/2)) - 3136.0)/(l_a^2*(-Yct*(-Yct^2/4 - 1.5) + Yct/2)*(-Yct*(-Yct^2/4 - 1.5)/2 + Yct/4)/(2*(Yct^2/4 + (-Yct^2/4 - 1.5)^2)) + 4.5*l_a^2*(-Yct*(-Yct^2/4 - 1.5)/2 + Yct/4)^2*sin(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)^2/(Yct^2/4 + (-Yct^2/4 - 1.5)^2) + 4.5*l_a^2*(-Yct*(-Yct^2/4 - 1.5)/2 + Yct/4)^2*cos(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)^2/(Yct^2/4 + (-Yct^2/4 - 1.5)^2) + 320))/(4.5*l_s^2*sin(thSt)^2 + 4.5*l_s^2*cos(thSt)^2 - 5.0625*(2*l_a*l_s*(-Yct*(-Yct^2/4 - 1.5)/2 + Yct/4)*sin(thSt)*cos(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) + 2*l_a*l_s*(-Yct*(-Yct^2/4 - 1.5)/2 + Yct/4)*sin(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)*cos(thSt)/sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2))^2/(l_a^2*(-Yct*(-Yct^2/4 - 1.5) + Yct/2)*(-Yct*(-Yct^2/4 - 1.5)/2 + Yct/4)/(2*(Yct^2/4 + (-Yct^2/4 - 1.5)^2)) + 4.5*l_a^2*(-Yct*(-Yct^2/4 - 1.5)/2 + Yct/4)^2*sin(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)^2/(Yct^2/4 + (-Yct^2/4 - 1.5)^2) + 4.5*l_a^2*(-Yct*(-Yct^2/4 - 1.5)/2 + Yct/4)^2*cos(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)^2/(Yct^2/4 + (-Yct^2/4 - 1.5)^2) + 320)));
end


function [Ia] = Ia(Yct, Ycdt, thSt, thSdt, l_a, l_s) 
	Ia = (l_a^2);
end


function [Xproj] = Xproj(Yct, Ycdt, thSt, thSdt, l_a, l_s) 
	Xproj = (-l_a*sin(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5) + l_s*cos(thSt));
end


function [Yproj] = Yproj(Yct, Ycdt, thSt, thSdt, l_a, l_s) 
	Yproj = (-l_a*cos(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5) + l_s*sin(thSt));
end


function [L] = L(Yct, Ycdt, thSt, thSdt, l_a, l_s) 
	L = (160*Ycdt^2 - 3136.0*Yct + 58.8*l_a^2*cos(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5) + l_a^2*(-Ycdt*Yct*(-Yct^2/4 - 1.5)/2 + Ycdt*Yct/4)^2/(2*(Yct^2/4 + (-Yct^2/4 - 1.5)^2)) + 44.1*l_a*cos(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5) - 44.1*l_s*sin(thSt) + 2.25*(l_a*(-Ycdt*Yct*(-Yct^2/4 - 1.5)/2 + Ycdt*Yct/4)*sin(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) + l_s*thSdt*cos(thSt))^2 + 2.25*(-l_a*(-Ycdt*Yct*(-Yct^2/4 - 1.5)/2 + Ycdt*Yct/4)*cos(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - l_s*thSdt*sin(thSt))^2 - 11477.76);
end


function [T] = T(Yct, Ycdt, thSt, thSdt, l_a, l_s) 
	T = (160*Ycdt^2 + l_a^2*(-Ycdt*Yct*(-Yct^2/4 - 1.5)/2 + Ycdt*Yct/4)^2/(2*(Yct^2/4 + (-Yct^2/4 - 1.5)^2)) + 2.25*(l_a*(-Ycdt*Yct*(-Yct^2/4 - 1.5)/2 + Ycdt*Yct/4)*sin(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) + l_s*thSdt*cos(thSt))^2 + 2.25*(-l_a*(-Ycdt*Yct*(-Yct^2/4 - 1.5)/2 + Ycdt*Yct/4)*cos(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - l_s*thSdt*sin(thSt))^2);
end


function [armTipPosX] = armTipPosX(Yct, Ycdt, thSt, thSdt, l_a, l_s) 
	armTipPosX = (-l_a*sin(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5));
end


function [projSpeedSq] = projSpeedSq(Yct, Ycdt, thSt, thSdt, l_a, l_s) 
	projSpeedSq = ((l_a*(-Ycdt*Yct*(-Yct^2/4 - 1.5)/2 + Ycdt*Yct/4)*sin(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) + l_s*thSdt*cos(thSt))^2 + (-l_a*(-Ycdt*Yct*(-Yct^2/4 - 1.5)/2 + Ycdt*Yct/4)*cos(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - l_s*thSdt*sin(thSt))^2);
end


function [Ma] = Ma(Yct, Ycdt, thSt, thSdt, l_a, l_s) 
	Ma = (12*l_a);
end


function [T_c] = T_c(Yct, Ycdt, thSt, thSdt, l_a, l_s) 
	T_c = (160*Ycdt^2);
end


function [projVelX] = projVelX(Yct, Ycdt, thSt, thSdt, l_a, l_s) 
	projVelX = (-l_a*(-Ycdt*Yct*(-Yct^2/4 - 1.5)/2 + Ycdt*Yct/4)*cos(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - l_s*thSdt*sin(thSt));
end


function [V_p] = V_p(Yct, Ycdt, thSt, thSdt, l_a, l_s) 
	V_p = (-44.1*l_a*cos(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5) + 44.1*l_s*sin(thSt));
end


function [projPos] = projPos(Yct, Ycdt, thSt, thSdt, l_a, l_s) 
	projPos = (MutableDenseMatrix([[-l_a*sin(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5) + l_s*cos(thSt)], [-l_a*cos(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5) + l_s*sin(thSt)]]));
end


function [Xcam] = Xcam(Yct, Ycdt, thSt, thSdt, l_a, l_s) 
	Xcam = (-3.33000000000000);
end


function [Xpb] = Xpb(Yct, Ycdt, thSt, thSdt, l_a, l_s) 
	Xpb = (Yct^2/4 - 1.83);
end


function [V_armcg] = V_armcg(Yct, Ycdt, thSt, thSdt, l_a, l_s) 
	V_armcg = (-58.8*l_a^2*cos(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5));
end


function [Ycdd] = Ycdd(Yct, Ycdt, thSt, thSdt, l_a, l_s) 
	Ycdd = ((-l_a^2*(-Yct*(-Yct^2/4 - 1.5) + Yct/2)*(Ycdt^2*Yct^2/4 - Ycdt^2*(-Yct^2/4 - 1.5)/2 + Ycdt^2/4)/(2*(Yct^2/4 + (-Yct^2/4 - 1.5)^2)) - l_a^2*(-Yct*(-Yct^2/4 - 1.5) + Yct/2)*(-Ycdt*Yct*(-Yct^2/4 - 1.5)/2 + Ycdt*Yct/4)*(Ycdt*Yct*(-Yct^2/4 - 1.5) - Ycdt*Yct/2)/(2*(Yct^2/4 + (-Yct^2/4 - 1.5)^2)^2) + l_a^2*(Yct*(-Yct^2/4 - 1.5) - Yct/2)*(-Ycdt*Yct*(-Yct^2/4 - 1.5)/2 + Ycdt*Yct/4)^2/(2*(Yct^2/4 + (-Yct^2/4 - 1.5)^2)^2) - 58.8*l_a^2*(-Yct*(-Yct^2/4 - 1.5)/2 + Yct/4)*sin(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 4.5*l_a*(-Yct*(-Yct^2/4 - 1.5)/2 + Yct/4)*(-Ycdt*Yct*(-Yct^2/4 - 1.5)/2 + Ycdt*Yct/4)*(l_a*(-Ycdt*Yct*(-Yct^2/4 - 1.5)/2 + Ycdt*Yct/4)*sin(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) + l_s*thSdt*cos(thSt))*cos(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 4.5*l_a*(-Yct*(-Yct^2/4 - 1.5)/2 + Yct/4)*(-Ycdt*Yct*(-Yct^2/4 - 1.5)/2 + Ycdt*Yct/4)*(-l_a*(-Ycdt*Yct*(-Yct^2/4 - 1.5)/2 + Ycdt*Yct/4)*cos(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - l_s*thSdt*sin(thSt))*sin(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/(Yct^2/4 + (-Yct^2/4 - 1.5)^2) + 4.5*l_a*(-Yct*(-Yct^2/4 - 1.5)/2 + Yct/4)*(l_a*(-Ycdt*Yct*(-Yct^2/4 - 1.5)/2 + Ycdt*Yct/4)^2*sin(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - l_a*(Ycdt^2*Yct^2/4 - Ycdt^2*(-Yct^2/4 - 1.5)/2 + Ycdt^2/4)*cos(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - l_a*(-Ycdt*Yct*(-Yct^2/4 - 1.5)/2 + Ycdt*Yct/4)*(Ycdt*Yct*(-Yct^2/4 - 1.5)/2 - Ycdt*Yct/4)*cos(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/(Yct^2/4 + (-Yct^2/4 - 1.5)^2)^(3/2) - l_s*thSdt^2*cos(thSt))*cos(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 4.5*l_a*(-Yct*(-Yct^2/4 - 1.5)/2 + Yct/4)*(l_a*(-Ycdt*Yct*(-Yct^2/4 - 1.5)/2 + Ycdt*Yct/4)^2*cos(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/(Yct^2/4 + (-Yct^2/4 - 1.5)^2) + l_a*(Ycdt^2*Yct^2/4 - Ycdt^2*(-Yct^2/4 - 1.5)/2 + Ycdt^2/4)*sin(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) + l_a*(-Ycdt*Yct*(-Yct^2/4 - 1.5)/2 + Ycdt*Yct/4)*(Ycdt*Yct*(-Yct^2/4 - 1.5)/2 - Ycdt*Yct/4)*sin(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/(Yct^2/4 + (-Yct^2/4 - 1.5)^2)^(3/2) - l_s*thSdt^2*sin(thSt))*sin(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 44.1*l_a*(-Yct*(-Yct^2/4 - 1.5)/2 + Yct/4)*sin(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 4.5*l_a*(l_a*(-Ycdt*Yct*(-Yct^2/4 - 1.5)/2 + Ycdt*Yct/4)*sin(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) + l_s*thSdt*cos(thSt))*(Ycdt*Yct^2/4 - Ycdt*(-Yct^2/4 - 1.5)/2 + Ycdt/4)*sin(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) + 4.5*l_a*(-l_a*(-Ycdt*Yct*(-Yct^2/4 - 1.5)/2 + Ycdt*Yct/4)*cos(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - l_s*thSdt*sin(thSt))*(Ycdt*Yct^2/4 - Ycdt*(-Yct^2/4 - 1.5)/2 + Ycdt/4)*cos(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 4.5*l_a*(-Yct*(-Yct^2/4 - 1.5)/2 + Yct/4)*(Ycdt*Yct*(-Yct^2/4 - 1.5)/2 - Ycdt*Yct/4)*(l_a*(-Ycdt*Yct*(-Yct^2/4 - 1.5)/2 + Ycdt*Yct/4)*sin(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) + l_s*thSdt*cos(thSt))*sin(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/(Yct^2/4 + (-Yct^2/4 - 1.5)^2)^(3/2) + 4.5*l_a*(-Yct*(-Yct^2/4 - 1.5)/2 + Yct/4)*(Ycdt*Yct*(-Yct^2/4 - 1.5)/2 - Ycdt*Yct/4)*(-l_a*(-Ycdt*Yct*(-Yct^2/4 - 1.5)/2 + Ycdt*Yct/4)*cos(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - l_s*thSdt*sin(thSt))*cos(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/(Yct^2/4 + (-Yct^2/4 - 1.5)^2)^(3/2) + 2.25*(l_a*(-Ycdt*Yct*(-Yct^2/4 - 1.5)/2 + Ycdt*Yct/4)*sin(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) + l_s*thSdt*cos(thSt))*(2*l_a*(-Yct*(-Yct^2/4 - 1.5)/2 + Yct/4)*(-Ycdt*Yct*(-Yct^2/4 - 1.5)/2 + Ycdt*Yct/4)*cos(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/(Yct^2/4 + (-Yct^2/4 - 1.5)^2) + 2*l_a*(Ycdt*Yct^2/4 - Ycdt*(-Yct^2/4 - 1.5)/2 + Ycdt/4)*sin(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) + 2*l_a*(Yct*(-Yct^2/4 - 1.5)/2 - Yct/4)*(-Ycdt*Yct*(-Yct^2/4 - 1.5)/2 + Ycdt*Yct/4)*sin(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/(Yct^2/4 + (-Yct^2/4 - 1.5)^2)^(3/2)) + 2.25*(-l_a*(-Ycdt*Yct*(-Yct^2/4 - 1.5)/2 + Ycdt*Yct/4)*cos(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - l_s*thSdt*sin(thSt))*(2*l_a*(-Yct*(-Yct^2/4 - 1.5)/2 + Yct/4)*(-Ycdt*Yct*(-Yct^2/4 - 1.5)/2 + Ycdt*Yct/4)*sin(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 2*l_a*(Ycdt*Yct^2/4 - Ycdt*(-Yct^2/4 - 1.5)/2 + Ycdt/4)*cos(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 2*l_a*(Yct*(-Yct^2/4 - 1.5)/2 - Yct/4)*(-Ycdt*Yct*(-Yct^2/4 - 1.5)/2 + Ycdt*Yct/4)*cos(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/(Yct^2/4 + (-Yct^2/4 - 1.5)^2)^(3/2)) - 2.25*(2*l_a*l_s*(-Yct*(-Yct^2/4 - 1.5)/2 + Yct/4)*sin(thSt)*cos(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) + 2*l_a*l_s*(-Yct*(-Yct^2/4 - 1.5)/2 + Yct/4)*sin(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)*cos(thSt)/sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2))*(4.5*l_s*(l_a*(-Ycdt*Yct*(-Yct^2/4 - 1.5)/2 + Ycdt*Yct/4)^2*sin(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - l_a*(Ycdt^2*Yct^2/4 - Ycdt^2*(-Yct^2/4 - 1.5)/2 + Ycdt^2/4)*cos(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - l_a*(-Ycdt*Yct*(-Yct^2/4 - 1.5)/2 + Ycdt*Yct/4)*(Ycdt*Yct*(-Yct^2/4 - 1.5)/2 - Ycdt*Yct/4)*cos(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/(Yct^2/4 + (-Yct^2/4 - 1.5)^2)^(3/2) - l_s*thSdt^2*cos(thSt))*sin(thSt) - 4.5*l_s*(l_a*(-Ycdt*Yct*(-Yct^2/4 - 1.5)/2 + Ycdt*Yct/4)^2*cos(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/(Yct^2/4 + (-Yct^2/4 - 1.5)^2) + l_a*(Ycdt^2*Yct^2/4 - Ycdt^2*(-Yct^2/4 - 1.5)/2 + Ycdt^2/4)*sin(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) + l_a*(-Ycdt*Yct*(-Yct^2/4 - 1.5)/2 + Ycdt*Yct/4)*(Ycdt*Yct*(-Yct^2/4 - 1.5)/2 - Ycdt*Yct/4)*sin(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/(Yct^2/4 + (-Yct^2/4 - 1.5)^2)^(3/2) - l_s*thSdt^2*sin(thSt))*cos(thSt) - 44.1*l_s*cos(thSt) - 2.25*(2*l_a*l_s*(-Yct*(-Yct^2/4 - 1.5)/2 + Yct/4)*sin(thSt)*cos(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) + 2*l_a*l_s*(-Yct*(-Yct^2/4 - 1.5)/2 + Yct/4)*sin(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)*cos(thSt)/sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2))*(-l_a^2*(-Yct*(-Yct^2/4 - 1.5) + Yct/2)*(Ycdt^2*Yct^2/4 - Ycdt^2*(-Yct^2/4 - 1.5)/2 + Ycdt^2/4)/(2*(Yct^2/4 + (-Yct^2/4 - 1.5)^2)) - l_a^2*(-Yct*(-Yct^2/4 - 1.5) + Yct/2)*(-Ycdt*Yct*(-Yct^2/4 - 1.5)/2 + Ycdt*Yct/4)*(Ycdt*Yct*(-Yct^2/4 - 1.5) - Ycdt*Yct/2)/(2*(Yct^2/4 + (-Yct^2/4 - 1.5)^2)^2) + l_a^2*(Yct*(-Yct^2/4 - 1.5) - Yct/2)*(-Ycdt*Yct*(-Yct^2/4 - 1.5)/2 + Ycdt*Yct/4)^2/(2*(Yct^2/4 + (-Yct^2/4 - 1.5)^2)^2) - 58.8*l_a^2*(-Yct*(-Yct^2/4 - 1.5)/2 + Yct/4)*sin(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 4.5*l_a*(-Yct*(-Yct^2/4 - 1.5)/2 + Yct/4)*(-Ycdt*Yct*(-Yct^2/4 - 1.5)/2 + Ycdt*Yct/4)*(l_a*(-Ycdt*Yct*(-Yct^2/4 - 1.5)/2 + Ycdt*Yct/4)*sin(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) + l_s*thSdt*cos(thSt))*cos(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 4.5*l_a*(-Yct*(-Yct^2/4 - 1.5)/2 + Yct/4)*(-Ycdt*Yct*(-Yct^2/4 - 1.5)/2 + Ycdt*Yct/4)*(-l_a*(-Ycdt*Yct*(-Yct^2/4 - 1.5)/2 + Ycdt*Yct/4)*cos(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - l_s*thSdt*sin(thSt))*sin(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/(Yct^2/4 + (-Yct^2/4 - 1.5)^2) + 4.5*l_a*(-Yct*(-Yct^2/4 - 1.5)/2 + Yct/4)*(l_a*(-Ycdt*Yct*(-Yct^2/4 - 1.5)/2 + Ycdt*Yct/4)^2*sin(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - l_a*(Ycdt^2*Yct^2/4 - Ycdt^2*(-Yct^2/4 - 1.5)/2 + Ycdt^2/4)*cos(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - l_a*(-Ycdt*Yct*(-Yct^2/4 - 1.5)/2 + Ycdt*Yct/4)*(Ycdt*Yct*(-Yct^2/4 - 1.5)/2 - Ycdt*Yct/4)*cos(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/(Yct^2/4 + (-Yct^2/4 - 1.5)^2)^(3/2) - l_s*thSdt^2*cos(thSt))*cos(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 4.5*l_a*(-Yct*(-Yct^2/4 - 1.5)/2 + Yct/4)*(l_a*(-Ycdt*Yct*(-Yct^2/4 - 1.5)/2 + Ycdt*Yct/4)^2*cos(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/(Yct^2/4 + (-Yct^2/4 - 1.5)^2) + l_a*(Ycdt^2*Yct^2/4 - Ycdt^2*(-Yct^2/4 - 1.5)/2 + Ycdt^2/4)*sin(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) + l_a*(-Ycdt*Yct*(-Yct^2/4 - 1.5)/2 + Ycdt*Yct/4)*(Ycdt*Yct*(-Yct^2/4 - 1.5)/2 - Ycdt*Yct/4)*sin(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/(Yct^2/4 + (-Yct^2/4 - 1.5)^2)^(3/2) - l_s*thSdt^2*sin(thSt))*sin(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 44.1*l_a*(-Yct*(-Yct^2/4 - 1.5)/2 + Yct/4)*sin(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 4.5*l_a*(l_a*(-Ycdt*Yct*(-Yct^2/4 - 1.5)/2 + Ycdt*Yct/4)*sin(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) + l_s*thSdt*cos(thSt))*(Ycdt*Yct^2/4 - Ycdt*(-Yct^2/4 - 1.5)/2 + Ycdt/4)*sin(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) + 4.5*l_a*(-l_a*(-Ycdt*Yct*(-Yct^2/4 - 1.5)/2 + Ycdt*Yct/4)*cos(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - l_s*thSdt*sin(thSt))*(Ycdt*Yct^2/4 - Ycdt*(-Yct^2/4 - 1.5)/2 + Ycdt/4)*cos(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 4.5*l_a*(-Yct*(-Yct^2/4 - 1.5)/2 + Yct/4)*(Ycdt*Yct*(-Yct^2/4 - 1.5)/2 - Ycdt*Yct/4)*(l_a*(-Ycdt*Yct*(-Yct^2/4 - 1.5)/2 + Ycdt*Yct/4)*sin(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) + l_s*thSdt*cos(thSt))*sin(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/(Yct^2/4 + (-Yct^2/4 - 1.5)^2)^(3/2) + 4.5*l_a*(-Yct*(-Yct^2/4 - 1.5)/2 + Yct/4)*(Ycdt*Yct*(-Yct^2/4 - 1.5)/2 - Ycdt*Yct/4)*(-l_a*(-Ycdt*Yct*(-Yct^2/4 - 1.5)/2 + Ycdt*Yct/4)*cos(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - l_s*thSdt*sin(thSt))*cos(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/(Yct^2/4 + (-Yct^2/4 - 1.5)^2)^(3/2) + 2.25*(l_a*(-Ycdt*Yct*(-Yct^2/4 - 1.5)/2 + Ycdt*Yct/4)*sin(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) + l_s*thSdt*cos(thSt))*(2*l_a*(-Yct*(-Yct^2/4 - 1.5)/2 + Yct/4)*(-Ycdt*Yct*(-Yct^2/4 - 1.5)/2 + Ycdt*Yct/4)*cos(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/(Yct^2/4 + (-Yct^2/4 - 1.5)^2) + 2*l_a*(Ycdt*Yct^2/4 - Ycdt*(-Yct^2/4 - 1.5)/2 + Ycdt/4)*sin(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) + 2*l_a*(Yct*(-Yct^2/4 - 1.5)/2 - Yct/4)*(-Ycdt*Yct*(-Yct^2/4 - 1.5)/2 + Ycdt*Yct/4)*sin(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/(Yct^2/4 + (-Yct^2/4 - 1.5)^2)^(3/2)) + 2.25*(-l_a*(-Ycdt*Yct*(-Yct^2/4 - 1.5)/2 + Ycdt*Yct/4)*cos(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - l_s*thSdt*sin(thSt))*(2*l_a*(-Yct*(-Yct^2/4 - 1.5)/2 + Yct/4)*(-Ycdt*Yct*(-Yct^2/4 - 1.5)/2 + Ycdt*Yct/4)*sin(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 2*l_a*(Ycdt*Yct^2/4 - Ycdt*(-Yct^2/4 - 1.5)/2 + Ycdt/4)*cos(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 2*l_a*(Yct*(-Yct^2/4 - 1.5)/2 - Yct/4)*(-Ycdt*Yct*(-Yct^2/4 - 1.5)/2 + Ycdt*Yct/4)*cos(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/(Yct^2/4 + (-Yct^2/4 - 1.5)^2)^(3/2)) - 3136.0)/(l_a^2*(-Yct*(-Yct^2/4 - 1.5) + Yct/2)*(-Yct*(-Yct^2/4 - 1.5)/2 + Yct/4)/(2*(Yct^2/4 + (-Yct^2/4 - 1.5)^2)) + 4.5*l_a^2*(-Yct*(-Yct^2/4 - 1.5)/2 + Yct/4)^2*sin(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)^2/(Yct^2/4 + (-Yct^2/4 - 1.5)^2) + 4.5*l_a^2*(-Yct*(-Yct^2/4 - 1.5)/2 + Yct/4)^2*cos(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)^2/(Yct^2/4 + (-Yct^2/4 - 1.5)^2) + 320))/(4.5*l_s^2*sin(thSt)^2 + 4.5*l_s^2*cos(thSt)^2 - 5.0625*(2*l_a*l_s*(-Yct*(-Yct^2/4 - 1.5)/2 + Yct/4)*sin(thSt)*cos(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) + 2*l_a*l_s*(-Yct*(-Yct^2/4 - 1.5)/2 + Yct/4)*sin(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)*cos(thSt)/sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2))^2/(l_a^2*(-Yct*(-Yct^2/4 - 1.5) + Yct/2)*(-Yct*(-Yct^2/4 - 1.5)/2 + Yct/4)/(2*(Yct^2/4 + (-Yct^2/4 - 1.5)^2)) + 4.5*l_a^2*(-Yct*(-Yct^2/4 - 1.5)/2 + Yct/4)^2*sin(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)^2/(Yct^2/4 + (-Yct^2/4 - 1.5)^2) + 4.5*l_a^2*(-Yct*(-Yct^2/4 - 1.5)/2 + Yct/4)^2*cos(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)^2/(Yct^2/4 + (-Yct^2/4 - 1.5)^2) + 320)) - 3136.0)/(l_a^2*(-Yct*(-Yct^2/4 - 1.5) + Yct/2)*(-Yct*(-Yct^2/4 - 1.5)/2 + Yct/4)/(2*(Yct^2/4 + (-Yct^2/4 - 1.5)^2)) + 4.5*l_a^2*(-Yct*(-Yct^2/4 - 1.5)/2 + Yct/4)^2*sin(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)^2/(Yct^2/4 + (-Yct^2/4 - 1.5)^2) + 4.5*l_a^2*(-Yct*(-Yct^2/4 - 1.5)/2 + Yct/4)^2*cos(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)^2/(Yct^2/4 + (-Yct^2/4 - 1.5)^2) + 320));
end


function [armTipPosY] = armTipPosY(Yct, Ycdt, thSt, thSdt, l_a, l_s) 
	armTipPosY = (-l_a*cos(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5));
end


function [armTipPos] = armTipPos(Yct, Ycdt, thSt, thSdt, l_a, l_s) 
	armTipPos = (MutableDenseMatrix([[-l_a*sin(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)], [-l_a*cos(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)]]));
end


function [projVelY] = projVelY(Yct, Ycdt, thSt, thSdt, l_a, l_s) 
	projVelY = (l_a*(-Ycdt*Yct*(-Yct^2/4 - 1.5)/2 + Ycdt*Yct/4)*sin(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5)/sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) + l_s*thSdt*cos(thSt));
end


function [T_a] = T_a(Yct, Ycdt, thSt, thSdt, l_a, l_s) 
	T_a = (l_a^2*(-Ycdt*Yct*(-Yct^2/4 - 1.5)/2 + Ycdt*Yct/4)^2/(2*(Yct^2/4 + (-Yct^2/4 - 1.5)^2)));
end


function [V_c] = V_c(Yct, Ycdt, thSt, thSdt, l_a, l_s) 
	V_c = (3136.0*Yct + 11477.76);
end


function [V] = V(Yct, Ycdt, thSt, thSdt, l_a, l_s) 
	V = (3136.0*Yct - 58.8*l_a^2*cos(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5) - 44.1*l_a*cos(sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - 1.5) + 44.1*l_s*sin(thSt) + 11477.76);
end


function [theta_arm] = theta_arm(Yct, Ycdt, thSt, thSdt, l_a, l_s) 
	theta_arm = (-sqrt(Yct^2/4 + (-Yct^2/4 - 1.5)^2) - pi/2 + 1.5);
end


function [Ypb] = Ypb(Yct, Ycdt, thSt, thSdt, l_a, l_s) 
	Ypb = (Yct/2);
end


end
end
