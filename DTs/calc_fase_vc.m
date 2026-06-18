function fase_vc = calc_fase_vc(Pvc, Hvc, Fluid, h_l_vc, h_g_vc)

import py.CoolProp.CoolProp.*

if Hvc <= h_l_vc
    fase_vc = 0;   % líquid

elseif Hvc >= h_g_vc
    fase_vc = 1;   % vapor
    
else
    fase_vc = 2;   % bifàsic
end

end