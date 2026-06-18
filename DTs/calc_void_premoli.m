function epsilon_g = calc_void_premoli(G, D, xg, rho_l, rho_g, mu_l, sigma)

% CALC_VOID_PREMOLI  
    % G      - velocitat màssica [kg/(m2·s)]
    % sigma  - tensió superficial [N/m]

    
    % CAS MONOFÀSIC amb protecció numèrica
    
    if xg <= 0
        epsilon_g = 1e-6;
    
    elseif xg >= 1
        epsilon_g = 1 - 1e-6;
    else
    
    % CAS BIFÀSIC
    g_acc = 9.81;  % acceleració gravitatòria [m/s2]

    % Reynolds del líquid
    Re_l = abs(G) * D / mu_l;

    % Weber del líquid
    We_l = G^2 * D / (sigma * rho_l * g_acc);

    % Coefs. F1 i F2 de Premoli
    F1 = 1.578 * Re_l^(-0.19) * (rho_l / rho_g)^0.22;
    F2 = 0.0273 * We_l * Re_l^(-0.51) * (rho_l / rho_g)^(-0.08);

    % Void fraction homogènia - Suposant velocitats iguals en ambdues components
    E_homog = 1 / (1 + (1 - xg) / xg * rho_g / rho_l);

    % Factor "y"
    y = E_homog / (1 - E_homog);

    % Relació de lliscament (Slip ratio)
    S = 1 + F1 * sqrt(abs(y / (1 + y * F2) - y * F2));

    % Void fraction final
    epsilon_g = 1 / (1 + (1 - xg) / xg * (rho_g / rho_l) * S);

    % Seguretat numèrica
    epsilon_g = max(epsilon_g, 1e-6);
    epsilon_g = min(epsilon_g, 1 - 1e-6);

    end

end
