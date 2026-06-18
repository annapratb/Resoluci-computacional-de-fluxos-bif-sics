function rho = calc_rho(P, h, Fluid, G, D, xg, rho_l, rho_g, mu_l, sigma)

import py.CoolProp.CoolProp.*

% Calcula la densitat segons la fase del fluid

% Saturació
h_L = PropsSI('H','P',P,'Q',0, Fluid);
h_V = PropsSI('H','P',P,'Q',1, Fluid);

rho_L = PropsSI('D','P',P,'Q',0, Fluid);
rho_V = PropsSI('D','P',P,'Q',1, Fluid);

%epsilon_g_in = 1e-6;
epsilon_g_in = calc_void_premoli(G, D, xg, rho_l, rho_g, mu_l, sigma);

% DETERMINACIÓ DE LA FASE
 
if h < h_L
    % Líquid subrefredat
    rho = PropsSI('D','P',P,'H',h, Fluid);

elseif h > h_V
    % Vapor sobreescalfat
    rho = PropsSI('D','P',P,'H',h, Fluid);

else 
    % Bifàsic
    rho = rho_V * epsilon_g_in + rho_L * (1 - epsilon_g_in);
end