function fase = calc_fase(P, H, Fluid, h_l, h_g)

import py.CoolProp.CoolProp.*

if H <= h_l
    fase = 0;   % líquid

elseif H >= h_g
    fase = 1;   % vapor
    
else
    fase = 2;   % bifàsic
end

end
