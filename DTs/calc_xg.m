function xg = calc_xg(P, h, Fluid, h_L, h_V)

% IMPORTANT: importar CoolProp dins la funció
import py.CoolProp.CoolProp.*

% DETERMINACIÓ DE LA FASE
if h <= h_L
    % Líquid (subrefredat o saturat líquid)
    xg = 0;
elseif h >= h_V
    % Vapor (sobreescalfat o saturat vapor)
    xg = 1;
else
    % Zona bifàsica → usar CoolProp
    xg = PropsSI('Q','P',P,'H',h,Fluid);

    % Seguretat numèrica (molt recomanable)
    xg = max(0, min(1, xg));
end
end