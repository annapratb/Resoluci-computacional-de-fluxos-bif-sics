function fase = calc_fase(P, H, Fluid)

import py.CoolProp.CoolProp.*

h_l = PropsSI('H','P',P,'Q',0,Fluid);
h_g = PropsSI('H','P',P,'Q',1,Fluid);

if H <= h_l
    fase = 0;   % líquid

elseif H >= h_g
    fase = 1;   % vapor
    
else
    fase = 2;   % bifàsic
end

end
