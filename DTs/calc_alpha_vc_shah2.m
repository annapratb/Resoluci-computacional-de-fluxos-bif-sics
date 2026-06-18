function alpha_vc = calc_alpha_vc_shah2(q, fase_vc, mdot_vc, mu_vc, D, lambda_vc, cp_vc, xg_vc, rho_l_vc, rho_g_vc, mu_l_vc, cp_l_vc, lambda_l_vc, hfg_vc)

if fase_vc == 0 || fase_vc == 1
    % MONOFÀSIC: Correlació de Gnielinski 
    alpha_vc = calc_alpha_gnielinski(mdot_vc, mu_vc, D, lambda_vc, cp_vc);
    
elseif fase_vc == 2
    % BIFÀSIC: Correlació de Shah
    
    % Constante de gravedad
    GRAV = 9.81;
    
    % Velocidad másica [kg/m2s]
    G = (4.0 * mdot_vc) / (pi * D^2);

    % Protección numérica para xg_vc (evitar divisiones por cero o valores negativos)
    if xg_vc <= 0
        xg_vc = 1e-6;
    elseif xg_vc >= 1
        xg_vc = 1 - 1e-6;
    end

    % 1. Cálculos de números adimensionales
    Rel = (abs(G) * (1.0 - xg_vc) * D) / mu_l_vc;
    Prl = (mu_l_vc * cp_l_vc) / lambda_l_vc;
    Frl = (abs(G) / rho_l_vc)^2 / (GRAV * D);
    
    % Nusselt para el líquido (Dittus-Boelter)
    Nu = 0.023 * (Rel^0.8) * (Prl^0.4);
    
    % Co: Número de Convección
    Co = ((1.0 - xg_vc) / xg_vc)^0.8 * sqrt(rho_g_vc / rho_l_vc);
    
    % Bo: Número de Ebullición
    Bo = q / (G * hfg_vc);
    if Bo <= 0
        Bo = 1e-10; % Protección si no hay flujo de calor
    end

    % 2. Coeficiente de transferencia de calor del líquido (Alfal)
    Alfal = (Nu * lambda_l_vc) / D;    

    % 3. Cálculo del valor de N
    if Frl >= 0.04
        N = Co;
    else
        N = 0.38 * (Frl^-0.3) * Co;
    end

    % 4. Cálculo del valor de F (Factor de corrección de Shah)
    if Bo >= 11.0e-4
        F = 14.7;
    else
        F = 15.43;
    end

    % 5. Cálculo del coeficiente Psi (Factor de incremento)
    % Psi_cb: Convective Boiling
    Psi_cb = 1.8 / (N^0.8);

    if N > 1.0
        % Región de Ebullición Nucleada
        if Bo > 0.3e-4
            Psi_nb = 230.0 * sqrt(Bo);  
        else
            Psi_nb = 1.0 + 46.0 * sqrt(Bo);
        end
        Psi = max(Psi_cb, Psi_nb);

    elseif N > 0.1 && N <= 1.0
        % Región de Burbujeo o Transición
        Psi_bs = F * sqrt(Bo) * exp(2.74 * (N^-0.1));
        Psi = max(Psi_cb, Psi_bs);
        
    else % N <= 0.1
        % Región de Convección Forzada Dominante
        Psi_bs = F * sqrt(Bo) * exp(2.47 * (N^-0.15));
        Psi = max(Psi_cb, Psi_bs);
    end

    % 6. Resultado final
    alpha_vc = Psi * Alfal;
    
end

end

%% FUNCIÓ: Correlació de Gnielinski

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
