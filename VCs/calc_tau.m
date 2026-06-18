function tau = calc_tau(fase_vc, fase, mdot, D, S, abs_rug, rho, mu, xg, rho_l, rho_g, mu_l, mu_g, sigma)

g = 9.81;

if fase_vc == 0 || fase_vc == 1 || fase == 0 || fase == 1
    %  CAS MONOFÀSIC  

    % Velocitat
    v  = mdot / (rho * S);

    % Nombre de Reynolds
    Re = rho * v * D / mu;

    % Factor de fricció de Fanning (empíric)
    if Re < 2000
        f = 16 / Re;
    elseif Re > 5e3 && Re < 3e4
        f = 0.079 / Re^0.25;
    else
        f = 0.096 / Re^0.25;
    end

    % Tensió de cisalla
    tau = f * rho * v^2 / 8;

else
    %  CAS BIFÀSIC 

%       fase_vc  – fase del VC (0=líquid, 1=vapor, 2=bifàsic)
%       mdot     – flux màssic [kg/s]
%       D        – diàmetre intern del tub [m]
%       S        – secció transversal [m²]
%       abs_rug  – rugositat absoluta [m]
%
%  monofàsic
%       rho      – densitat de la fase [kg/m³]       (rho_tp_vc)
%       mu       – viscositat dinàmica [Pa·s]         (mu_vc)
%
%   bifàsic
%       xg       – qualitat (títol de vapor) [-]      (xg_vc)
%       rho_l    – densitat líquid saturat [kg/m³]    (rho_l_vc)
%       rho_g    – densitat vapor saturat [kg/m³]     (rho_g_vc)
%       mu_l     – viscositat dinàmica líquid [Pa·s]  (mu_l_vc)
%       mu_g     – viscositat dinàmica vapor [Pa·s]   (mu_g_vc)
%       sigma    – tensió superficial [N/m]           (sigma_vc)
   
    % Marcar límits del títol: si tenim el resultat anterior xg=-0.2 tria 0 i si tenim xig=1.05 tria 1.
    Xg = max(0, min(1, xg));


    % Avís de rang recomanat
    if mu_l / mu_g > 1000
        warning('calc_tau (Friedel): mu_l/mu_g > 1000 — fora del rang recomanat.');
    end

    % Flux màssic per unitat d'àrea  G = rho·v  [kg/(m²·s)]
    G = mdot / S;

    % Densitat homogènia RoH  (definició Friedel)
    RoH = (rho_g * rho_l) / (Xg * rho_l + (1 - Xg) * rho_g);

    % Densitat homogènia de la mescla, RoM
   % RoM=RoH;

    % Velocitat de la mescla  Vm
    %Vm = G / RoM;

    % Factors de fricció de Churchill (Darcy) per a cada fase
    fl = friction_Churchill(G, mu_l, D, abs_rug);
    fg = friction_Churchill(G, mu_g, D, abs_rug);

    % Nombres adimensionals de Weber i Froude
    We = G^2 * D / (RoH * sigma);
    Fr = G^2   / (g * D * RoH^2);

    % Coeficients E, F, H de Friedel
    E_f = (1 - Xg)^2 + Xg^2 * (rho_l * fg) / (rho_g * fl);
    F_f = Xg^0.78 * (1 - Xg)^0.224;
    H_f = (rho_l/rho_g)^0.91 * (mu_g/mu_l)^0.19 * (1 - mu_g/mu_l)^0.70;

    % Factor de correció de la caiguda de pressió
    Fi = sqrt(E_f + (3.23 * F_f * H_f) / (Fr^0.045 * We^0.035));

    % Tensió de cisalla bifàsica
    tau = (fl/4) * (G^2 / (2 * rho_l)) * Fi^2;

    % disp(['Xg = ', num2str(Xg), ' ']);
    % disp(['rho_g = ', num2str(rho_g), ' ']);
    % disp(['rho_l = ', num2str(rho_l), ' ']);
    % disp(['mu_g = ', num2str(mu_g), ' ']);
    % disp(['mu_l = ', num2str(mu_l), ' ']);
    % disp(['G = ', num2str(G), ' ']);
    % disp(['RoH = ', num2str(RoH), ' ']);
    % disp(['tau = ', num2str(tau), ' ']);


    % Si el valor de tau es indefinit (NaN) entra al if i para el programa
    if isnan(tau)
        error('calc_tau (Friedel): tau = NaN  (fl=%.4e, G=%.4e, rho_l=%.4e, Fi=%.4e)', ...
              fl, G, rho_l, Fi);
    end

    % Tensió de cisalla bifàsica
    %tau = 1;
end

end  % fi de calc_tau

%  SUBFUNCIÓ 

%    RoM     – densitat de la fase [kg/m³]
%    Vm      – velocitat de la fase [m/s]
%    visc    – viscositat dinàmica de la fase [Pa·s]
%    Ds      – diàmetre interior [m]
%    abs_rug – rugositat absoluta [m]

function f = friction_Churchill(G, visc, Ds, abs_rug)

Re = (G * Ds) / visc;

if Re == 0
    f = 0;
    return;
end

C_ch = (7 / Re)^0.9 + 0.27 * (abs_rug / Ds);
B_ch = (37530 / Re)^16;
A_ch = (2.457 * log(1 / C_ch))^16;   % log() = ln() en MATLAB

f = 8 * ((8 / Re)^12 + 1 / (A_ch + B_ch)^1.5)^(1/12);

if isnan(f)
    error('friction_Churchill: f = NaN  (Re=%.4e, A=%.4e, B=%.4e, C=%.4e)', ...
          Re, A_ch, B_ch, C_ch);
end

end  % fi de friction_Churchill