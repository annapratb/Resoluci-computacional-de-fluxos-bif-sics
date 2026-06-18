function rho = calc_rho(P, h, Fluid, h_L, h_V, rho_L_sat, rho_V_sat, eps)

import py.CoolProp.CoolProp.*

% Calcula la densitat segons la fase del fluid

% DETERMINACIÓ DE LA FASE
 
if h < h_L
    % Líquid subrefredat
    rho = PropsSI('D','P',P,'H',h, Fluid);

elseif h > h_V
    % Vapor sobreescalfat
    rho = PropsSI('D','P',P,'H',h, Fluid);

else 
    % Bifàsic
    rho = rho_V_sat * eps + rho_L_sat * (1 - eps);
end