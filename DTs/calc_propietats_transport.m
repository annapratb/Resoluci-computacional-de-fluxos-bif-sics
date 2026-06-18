function [mu, lambda, cp] = calc_propietats_transport(P, H, Fluid)
    % Evitar error CoolProp

    import py.CoolProp.CoolProp.PropsSI
    h_L = PropsSI('H','P',P,'Q',0, Fluid);
    h_V = PropsSI('H','P',P,'Q',1, Fluid);
    xg = calc_xg(P, H, Fluid, h_L, h_V);

    try
        if xg <= 0  % Líquid subrefredat
            mu = double(PropsSI('V', 'P', P, 'H', H, Fluid));
            lambda = double(PropsSI('L', 'P', P, 'H', H, Fluid));
            cp = double(PropsSI('C','P', P, 'H', H,Fluid));
            
        elseif xg >= 1 % Vapor sobreescalfat
            mu = double(PropsSI('V', 'P', P, 'H', H, Fluid));
            lambda = double(PropsSI('L', 'P', P, 'H', H, Fluid));
            cp = double(PropsSI('C','P', P, 'H', H, Fluid));
            
        else % ZONA BIFÀSICA (Interpolació lineal)
            mu_l = double(PropsSI('V', 'P', P, 'Q', 0, Fluid));
            mu_g = double(PropsSI('V', 'P', P, 'Q', 1, Fluid));
            mu = mu_l * (1 - xg) + mu_g * xg;
            
            lambda_l = double(PropsSI('L', 'P', P, 'Q', 0, Fluid));
            lambda_g = double(PropsSI('L', 'P', P, 'Q', 1, Fluid));
            lambda = lambda_l * (1 - xg) + lambda_g * xg;

            cp_l = double(PropsSI('C','P', P,'Q', 0, Fluid));
            cp_g = double(PropsSI('C','P', P,'Q', 1, Fluid));
            cp = cp_l * (1 - xg) + cp_g * xg;

            
        end
    catch
        % Si CoolProp falla fins i tot en zones monofàsiques prop de saturació
        % retornem propietats de saturació 
        if xg < 0.5
            mu = double(PropsSI('V', 'P', P, 'Q', 0, Fluid));
            lambda = double(PropsSI('L', 'P', P, 'Q', 0, Fluid));
            cp = double(PropsSI('C','P', P,'Q', 0, Fluid));
        else
            mu = double(PropsSI('V', 'P', P, 'Q', 1, Fluid));
            lambda = double(PropsSI('L', 'P', P, 'Q', 1, Fluid));
            cp = double(PropsSI('C','P', P,'Q', 1, Fluid));
        end
    end
end
