function alpha_vc = calc_alpha_vc(q, fase_vc, mdot_vc, mu_vc, D, lambda_vc, cp_vc, xg_vc, rho_l_vc, rho_g_vc, mu_l_vc, cp_l_vc, lambda_l_vc, hfg_vc)

if fase_vc == 0 || fase_vc == 1
    % MONOFÀSIC: Correlació de Gnielinski 
    
    alpha_vc = calc_alpha_gnielinski(mdot_vc, mu_vc, D, lambda_vc, cp_vc);
    
elseif fase_vc == 2
    % BIFÀSIC: Correlació de Gungor-Winterton
    
    % Velocitat màssica
    G = (4.0 * mdot_vc) / (pi * D^2);

    % Boiling number
    Bo = q / (G * hfg_vc);
    if Bo <= 0
        Bo = 1e-10; % Protecció si no hi ha flux de calor positiu (evita complexos)
    end

    % Froude
    Fr = G^2 / (rho_l_vc^2 * 9.81 * D);

    % Protecció numèrica: evitar divisió per zero quan xg_vc -> 1
    if (1 - xg_vc) < 1e-10
        xg_vc = 1 - 1e-10;
    end

    % Reynolds líquid
    Re_l = G * (1 - xg_vc) * D / mu_l_vc;

    % Prandtl líquid
    Pr_l = cp_l_vc * mu_l_vc / lambda_l_vc;

    % Factor de millora
    E = 1 + 3000 * Bo^0.86 + 1.12 * (xg_vc/(1-xg_vc))^0.75 * (rho_l_vc/rho_g_vc)^0.41;

    % Correcció per Froude baix
    if Fr < 0.05
        E = E * Fr^(0.1 - 2*Fr);
    end

    % Coeficient convectiu del líquid
    hl = 0.023 * Re_l^0.8 * Pr_l^0.4 * (lambda_l_vc / D);

    % Resultat final
    alpha_vc = E * hl;
    
end

end


%% FUNCIÓ: Correlació de Gnielinski
% Calcula el coeficient convectiu monofàsic a partir de:
%   mdot   - cabal màssic [kg/s]
%   mu     - viscositat dinàmica [Pa·s]
%   D      - diàmetre del tub [m]
%   lambda - conductivitat tèrmica [W/(m·K)]
%   cp     - calor específic [J/(kg·K)]

function alphaS = calc_alpha_gnielinski(mdot, mu, D, lambda, cp)

    % Velocitat màssica
    G = (4.0 * mdot) / (pi * D^2);
    
    % Reynolds
    Re = abs(G) * D / mu;
    
    if Re <= 2000
        Nu = 3.66;
    else
        % Factor de fricció 
        f = (1.82 * log10(Re) - 1.64)^(-2);
        
        % Prandtl
        Pr = (cp * mu) / lambda;
        
        % Nusselt 
        Nu = ((f/8)*(Re-1000)*Pr) / (1 + 12.7*sqrt(f/8)*(Pr^(2/3)-1));
    end
    
    alphaS = Nu * lambda / D;

end
