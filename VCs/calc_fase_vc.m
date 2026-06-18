function fase_vc = calc_fase_vc(Pvc, Hvc, Fluid)

import py.CoolProp.CoolProp.*

h_l_vc = PropsSI('H','P',Pvc,'Q',0,Fluid);
h_g_vc = PropsSI('H','P',Pvc,'Q',1,Fluid);

if Hvc <= h_l_vc
    fase_vc = 0;   % líquid

elseif Hvc >= h_g_vc
    fase_vc = 1;   % vapor
    
else
    fase_vc = 2;   % bifàsic
end

end