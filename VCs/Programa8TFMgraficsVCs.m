%% PROGRAMA 8 - Condicions (Estudi d'Independència de Malla)
clear all
close all
clc
import py.CoolProp.CoolProp.*

%% DADES INICIALITZACIÓ
Fluid = 'R134a';

%% GEOMETRIA I SUPOSICIONS GENERALS
L = 17.25;             % longitud total del tub [m]
D = 0.012;             % diàmetre interior [m]  (Kattan et al.)
S = pi * D^2/4;        % secció transversal [m²]
Perim = pi * D;        % perímetre interior [m]
abs_rug = 15e-7;       % rugositat absoluta del tub [m]
dt = 1e9;              % Temps molt gran (quasi estacionari)
g = 9.81;
theta = 0;             % tub horitzontal

Qdot = zeros(1, NUM);   
error_tol = 1e-5;
N = 100;


%% CONDICIÓ FRONTERA
tipus_condicio = 1;    % 1 --> flux de calor constant, 2 --> Text fixa

%% DEFINICIÓ DE LES MALLES A ESTUDIAR
llista_NUM = [10:10:50, 100:50:300];
n_malles = length(llista_NUM);

% Vectors per guardar els valors finals a la SORTIDA DEL TUB (node final) per a cada malla
P_sortida_malla = zeros(1, n_malles);
H_sortida_malla = zeros(1, n_malles);
T_sortida_malla = zeros(1, n_malles);
alpha_mitja_malla = zeros(1, n_malles); % Coeficient de transferència mitjà al tub


% Iteració sobre els diferents volums de control
for m = 1:n_malles
    NUM = llista_NUM(m);
    dz = L / NUM;         % longitud volum de control [m]
    A = Perim*dz;         % àrea mullada del VC [m²]

    fprintf('  EXECUCCIÓ AMB NUM = %d VOLUMS DE CONTROL\n', NUM);

    %% CONFIGURACIÓ DE FRONTERES SEGONS LA MALLA ACTUAL
    if tipus_condicio == 2 
        Text_ref = 210; 
        Text = Text_ref * ones(1, NUM);  
        P_in = 10^6;                     
        h_in = 1100000;                  
        mdot_in = 0.005;                 
    else
        q = 10000;                       
        Text = zeros(1, NUM);            
        T_sat_in = 10 + 273.15;           
        P_in = double(PropsSI('P','T',T_sat_in,'Q',0,Fluid)); 
        G = 300;                         
        mdot_in = G * S;                 
        h_l_in = double(PropsSI('H','P',P_in,'Q',0,Fluid));
        h_g_in = double(PropsSI('H','P',P_in,'Q',1,Fluid));
        h_in   = h_l_in + 1000;
    end

    %% Inicialització de Vectors de Propietats (mida dinàmica segons NUM)
    P = zeros(1, NUM+1);
    H = zeros(1, NUM+1);
    T = zeros(1, NUM+1);
    mdot = zeros(1, NUM+1);
    fase    = zeros(1, NUM+1);
    fase_vc = zeros(1, NUM);
    alpha_vc    = zeros(1, NUM);
    xg_vc_vec   = zeros(1, NUM);   
    
    P(1) = P_in; H(1) = h_in; mdot(1) = mdot_in;
    P(2:NUM+1) = P(1); H(2:NUM+1) = H(1); mdot(2:NUM+1) = mdot(1);

    rho_l = zeros(1, NUM+1); rho_g = zeros(1, NUM+1); cp = zeros(1, NUM+1);
    mu = zeros(1, NUM+1); lambda = zeros(1, NUM+1); epsilon_g = zeros(1, NUM+1);
    xg = zeros(1, NUM+1); rho_tp = zeros(1, NUM+1); v = zeros(1, NUM+1);

    for k = 1:NUM+1
        rho_l(k) = PropsSI('D','P',P(k),'Q',0,Fluid);
        rho_g(k) = PropsSI('D','P',P(k),'Q',1,Fluid);
        cp(k) = PropsSI('C','P',P(k),'H',H(k),Fluid);
        T(k) = PropsSI('T','P',P(k),'H',H(k),Fluid) - 273.15;
        [mu(k), lambda(k)] = calc_propietats_transport(P(k), H(k), Fluid); 
        xg(k) = calc_xg(P(k), H(k), Fluid);
        fase(k) = calc_fase(P(k), H(k), Fluid);  
        
        mu_l_k  = PropsSI('V','P',P(k),'Q',0,Fluid);
        sigma_k = PropsSI('I','P',P(k),'Q',0,Fluid);
        G_k = mdot(k) / S;
        epsilon_g(k) = calc_void_premoli(G_k, D, xg(k), rho_l(k), rho_g(k), mu_l_k, sigma_k);
        rho_tp(k) = calc_rho(P(k), H(k), Fluid, G_k, D, xg(k), rho_l(k), rho_g(k), mu_l_k, sigma_k);
        v(k) = mdot(k) / (rho_tp(k) * S);
    end

    %% DADES OLD (mida dinàmica segons NUM)
    rho_tp_old = rho_tp; xg_old = xg; mdot_old = mdot; rho_g_old = rho_g;
    rho_l_old = rho_l; epsilon_g_old = epsilon_g; P_old = P; H_old = H;

    rho_tp_vc_old = zeros(1, NUM); xg_old_vc = zeros(1, NUM); mdot_old_vc = zeros(1, NUM);
    rho_g_old_vc = zeros(1, NUM); rho_l_old_vc = zeros(1, NUM); epsilon_g_old_vc = zeros(1, NUM);
    P_old_vc = zeros(1, NUM); h_old_vc = zeros(1, NUM);

    for j = 1:NUM
        P_old_vc(j) = (P(j) + P(j+1)) / 2;
        h_old_vc(j) = (H(j) + H(j+1)) / 2;
        rho_l_old_vc(j) = PropsSI('D','P',P_old_vc(j),'Q',0,Fluid);
        rho_g_old_vc(j) = PropsSI('D','P',P_old_vc(j),'Q',1,Fluid);
        xg_old_vc(j) = calc_xg(P_old_vc(j), h_old_vc(j), Fluid);
        mdot_old_vc(j) = (mdot(j) + mdot(j+1)) / 2;
        
        mu_l_old_vc_j  = PropsSI('V','P',P_old_vc(j),'Q',0,Fluid);
        sigma_old_vc_j = PropsSI('I','P',P_old_vc(j),'Q',0,Fluid);
        G_vc_j = mdot_old_vc(j) / S;
        epsilon_g_old_vc(j) = calc_void_premoli(G_vc_j, D, xg_old_vc(j), rho_l_old_vc(j), rho_g_old_vc(j), mu_l_old_vc_j, sigma_old_vc_j);
        rho_tp_vc_old(j) = calc_rho(P_old_vc(j), h_old_vc(j), Fluid, G_vc_j, D, xg_old_vc(j), rho_l_old_vc(j), rho_g_old_vc(j), mu_l_old_vc_j, sigma_old_vc_j);
    end

    %% BUCLE ITERATIU (Resolució espacial per a la malla actual)
    for j = 1:NUM      
        P(j+1)=P(j); H(j+1)=H(j); xg(j+1)=xg(j);
        rho_l(j+1)=rho_l(j); rho_g(j+1)=rho_g(j);
        epsilon_g(j+1)=epsilon_g(j); alpha_vc(j)=1000;
        
        for i = 1:N    
            mdot_vc = (mdot(j) + mdot(j+1)) / 2;
            Pvc = (P(j) + P(j+1)) / 2;
            Hvc = (H(j) + H(j+1)) / 2;
            
            rho_l_vc = PropsSI('D','P',Pvc,'Q',0,Fluid);
            rho_g_vc = PropsSI('D','P',Pvc,'Q',1,Fluid);
            cp_vc = PropsSI('C','P',Pvc,'H',Hvc,Fluid);   
            [mu_vc, lambda_vc] = calc_propietats_transport(Pvc, Hvc, Fluid);
            Tvc = PropsSI('T','P',Pvc,'H',Hvc,Fluid) - 273.15;
            xg_vc  = calc_xg(Pvc, Hvc, Fluid);
            fase_vc(j) = calc_fase_vc(Pvc, Hvc, Fluid);  
            xg_vc_vec(j) = max(0, min(1, xg_vc));        
            
            mu_l_vc     = PropsSI('V','P',Pvc,'Q',0,Fluid);
            mu_g_vc     = PropsSI('V','P',Pvc,'Q',1,Fluid);    
            cp_l_vc     = PropsSI('C','P',Pvc,'Q',0,Fluid);
            lambda_l_vc = PropsSI('L','P',Pvc,'Q',0,Fluid);
            h_l_vc      = PropsSI('H','P',Pvc,'Q',0,Fluid);
            h_g_vc      = PropsSI('H','P',Pvc,'Q',1,Fluid);
            hfg_vc      = h_g_vc - h_l_vc;                     
            sigma_vc    = PropsSI('I','P',Pvc,'Q',0,Fluid);
            
            G_vc = mdot_vc / S;
            epsilon_g_vc = calc_void_premoli(G_vc, D, xg_vc, rho_l_vc, rho_g_vc, mu_l_vc, sigma_vc);
            rho_tp_vc = calc_rho(Pvc, Hvc, Fluid, G_vc, D, xg_vc, rho_l_vc, rho_g_vc, mu_l_vc, sigma_vc);
            v_vc = mdot_vc / (rho_tp_vc * S);
            
            alpha_vc(j) = calc_alpha_vc(Qdot(j)/A, fase_vc(j), mdot_vc, mu_vc, D, lambda_vc, cp_vc, xg_vc, rho_l_vc, rho_g_vc, mu_l_vc, cp_l_vc, lambda_l_vc, hfg_vc);
            
            P_ite = P(j+1); H_ite = H(j+1); mdot_ite = mdot(j+1);
            
            %% EQUACIÓ DE LA MASSA
            mdot(j+1) = mdot(j) - ((rho_tp(j) - rho_tp_old(j))/dt) * S * dz;
            
            %% EQUACIÓ MOMENTUM
            tau = calc_tau(fase_vc(j), fase(j), mdot_vc, D, S, abs_rug, rho_tp_vc, mu_vc, xg_vc, rho_l_vc, rho_g_vc, mu_l_vc, mu_g_vc, sigma_vc);
            term_1 = (mdot(j+1)/S)^2 * ( xg(j+1)^2/(rho_g(j+1)*epsilon_g(j+1)) + (1-xg(j+1))^2/(rho_l(j+1)*(1-epsilon_g(j+1))) );
            term_2 = (mdot(j)/S)^2 * ( xg(j)^2/(rho_g(j)*epsilon_g(j)) + (1-xg(j))^2/(rho_l(j)*(1-epsilon_g(j))) );
            term_t = (mdot_vc - mdot_old_vc(j))/(S*dt);
            term_3 = term_t + (tau * Perim / S) + rho_tp_vc * g * sin(theta);
            
            P(j+1) = P(j) - term_1 + term_2 - term_3 * dz;
            
            %% EQUACIÓ ENERGIA
            a = ( xg(j+1) / 2 ) * ( (mdot(j+1)^2 * xg(j+1)^2) / (rho_g(j+1)^2 * epsilon_g(j+1)^2 * S^2) - (mdot(j)^2 * xg(j)^2) / (rho_g(j)^2 * epsilon_g(j)^2 * S^2) ) + ...
                (( 1 - xg(j+1)) / 2) * ( (mdot(j+1)^2 * (1 - xg(j+1))^2) / (rho_l(j+1)^2 * (1 - epsilon_g(j+1))^2 * S^2) - (mdot(j)^2 * (1 - xg(j))^2) / (rho_l(j)^2 * ( 1 - epsilon_g(j) )^2 * S^2) ) ...
                + (g*sin(theta)*dz - H(j)); 
                
            b = ( xg(j) / 2) * ( (mdot(j)^2 * xg(j)^2) / (rho_g(j)^2 * epsilon_g(j)^2 * S^2) - (mdot(j+1)^2 * xg(j+1)^2) / (rho_g(j+1)^2 * epsilon_g(j+1)^2 * S^2) ) + ...
                ((1 - xg(j)) / 2) * ( (mdot(j)^2 * (1 - xg(j))^2) / (rho_l(j)^2 * (1 - epsilon_g(j))^2 * S^2) - (mdot(j+1)^2 * (1 - xg(j+1))^2) / (rho_l(j+1)^2 * (1 - epsilon_g(j+1) )^2 * S^2) ) ...
                + (-g*sin(theta)*dz + H(j));
                
            c = rho_tp_vc_old(j) * (H_old(j+1) + H_old(j) - H(j)) + ( (P(j) + P(j+1)) - (P_old(j) + P_old(j+1)) ) ...
                - (rho_g_old_vc(j) * epsilon_g_old_vc(j))/2 * ( (mdot(j)^2  * xg(j)^2) / (rho_g(j)^2 * epsilon_g(j)^2 * S^2) + (mdot(j+1)^2 * xg(j+1)^2) / (rho_g(j+1)^2 * epsilon_g(j+1)^2 * S^2) ) ...
                - (rho_l_old_vc(j) * (1 - epsilon_g_old_vc(j)))/2 * ( (mdot(j)^2  * (1 - xg(j))^2) / (rho_l(j)^2 * (1 - epsilon_g(j))^2 * S^2) + (mdot(j+1)^2 * (1 - xg(j+1))^2) / (rho_l(j+1)^2 * (1 - epsilon_g(j+1))^2 * S^2) ) ...
                + (rho_g_old_vc(j) * epsilon_g_old_vc(j))/2 * ( (mdot_old(j)^2 * xg_old(j)^2) / (rho_g_old(j)^2 * epsilon_g(j)^2 * S^2) + (mdot_old(j+1)^2 * xg_old(j+1)^2) / (rho_g_old(j+1)^2 * epsilon_g_old(j+1)^2 * S^2) ) ...
                + (rho_l_old_vc(j) * (1 - epsilon_g_old_vc(j)))/2 * ( (mdot_old(j)^2 * (1 - xg_old(j))^2) / (rho_l_old(j)^2 * (1 - epsilon_g_old(j))^2 * S^2) + (mdot_old(j+1)^2 * (1 - xg_old(j+1))^2) / (rho_l_old(j+1)^2 * (1 - epsilon_g_old(j+1))^2 * S^2) );
            
            if tipus_condicio == 1
                Qdot(j) = q * A;
                if alpha_vc(j) > 0
                    Text(j) = Qdot(j) / (alpha_vc(j) * A) + Tvc;
                else
                    Text(j) = Tvc;
                end
            else
                Qdot(j) = alpha_vc(j) * A * (Text(j) - Tvc);
            end
            
            H(j+1) = (2*Qdot(j) - mdot(j+1)*a + mdot(j)*b + ((S*dz/dt)*c)) / ((rho_tp_vc_old(j)*S*dz)/dt + mdot(j+1) + mdot(j));
            
            %% RECÀLCUL DE PROPIETATS Node (j+1)
            rho_l(j+1) = PropsSI('D','P',P(j+1),'Q',0, Fluid);
            rho_g(j+1) = PropsSI('D','P',P(j+1),'Q',1, Fluid);     
            cp(j+1) = PropsSI('C','P',P(j+1),'H',H(j+1), Fluid);
            T(j+1) = PropsSI('T','P',P(j+1),'H',H(j+1), Fluid) - 273.15;
            [mu(j+1), lambda(j+1)] = calc_propietats_transport(P(j+1), H(j+1), Fluid);
            xg(j+1) = calc_xg(P(j+1), H(j+1), Fluid); 
            fase(j+1) = calc_fase(P(j+1), H(j+1), Fluid);  
            
            mu_l_jp1  = PropsSI('V','P',P(j+1),'Q',0,Fluid);
            sigma_jp1 = PropsSI('I','P',P(j+1),'Q',0,Fluid);
            G_jp1 = mdot(j+1) / S;
            epsilon_g(j+1) = calc_void_premoli(G_jp1, D, xg(j+1), rho_l(j+1), rho_g(j+1), mu_l_jp1, sigma_jp1);
            rho_tp(j+1) = calc_rho(P(j+1), H(j+1), Fluid, G_jp1, D, xg(j+1), rho_l(j+1), rho_g(j+1), mu_l_jp1, sigma_jp1);
            v(j+1) = mdot(j+1) / (rho_tp(j+1) * S);
            
            e1 = abs(mdot(j+1) - mdot_ite);
            e2 = abs(P(j+1) - P_ite);
            e3 = abs(H(j+1) - H_ite);
            
            if (e1 < error_tol && e2 < error_tol && e3 < error_tol)
                break;
            end
        end % Fi bucle intern
  
    end % Fi bucle espacial (j)

    %% DESAR VALORS DE SORTIDA PER A LA MALLA ACTUAL
    % El node final és el NUM+1, corresponent a la sortida física del tub (z = L)
    P_sortida_malla(m) = P(NUM+1);
    H_sortida_malla(m) = H(NUM+1);
    T_sortida_malla(m) = T(NUM+1);
    alpha_mitja_malla(m) = mean(alpha_vc); % Fem la mitjana del coeficient de transferència al tub
    
end % Fi bucle malles (m)

%% MOSTRAR TAULA COMPARTIVA I GRÀFICS D'INDEPENDÈNCIA DE MALLA
disp(' ');
disp('============= RESUM DE CONVERGÈNCIA ESPACIAL (INDEPENDÈNCIA DE MALLA) =============');
fprintf('%6s | %14s | %14s | %12s | %16s\n', 'NUM VC', 'P Sortida (Pa)', 'H Sortida (J/kg)', 'T Sortida (ºC)', 'alpha mitjà (W/m²K)');
fprintf('----------------------------------------------------------------------------------------\n');
for m = 1:n_malles
    fprintf('%6d | %14.2f | %14.1f | %12.2f | %16.2f\n', ...
        llista_NUM(m), P_sortida_malla(m), H_sortida_malla(m), T_sortida_malla(m), alpha_mitja_malla(m));
end

% CODI MODIFICAT PER A FONS BLANC
% =========================================================================

% 1. Creem la figura i forcem que el seu fons general sigui blanc ('w')
fig = figure('Name', 'Estudi de Convergència Espacial', 'NumberTitle', 'off', 'Color', 'w');

% --- Subgràfic 1: Pressió ---
ax1 = subplot(2,2,1);
plot(llista_NUM, P_sortida_malla, 'k-s', 'LineWidth', 1.5, 'MarkerFaceColor', 'K');
grid on; 
xlabel('Número de Volums de Control'); ylabel('Pressió de Sortida [Pa]');
title('Evolució Pressió', 'Color', 'k');
set(ax1, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'GridColor', [0.5 0.5 0.5]);

% --- Subgràfic 2: Entalpia ---
ax2 = subplot(2,2,2);
plot(llista_NUM, H_sortida_malla, 'k-^', 'LineWidth', 1.5, 'MarkerFaceColor', 'k');
grid on; 
xlabel('Número de Volums de Control'); ylabel('Entalpia de Sortida [J/kg]');
title('Evolució Entalpia', 'Color', 'k');
set(ax2, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'GridColor', [0.5 0.5 0.5]);

% --- Subgràfic 3: Temperatura ---
ax3 = subplot(2,2,3);
plot(llista_NUM, T_sortida_malla, 'k-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'k');
grid on; 
xlabel('Número de Volums de Control'); ylabel('Temperatura de Sortida [ºC]');
title('Evolució Temperatura', 'Color', 'k');
set(ax3, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'GridColor', [0.5 0.5 0.5]);

% --- Subgràfic 4: Coef. Transferència ---
ax4 = subplot(2,2,4);
plot(llista_NUM, alpha_mitja_malla, 'k-d', 'LineWidth', 1.5, 'MarkerFaceColor', 'k');
grid on; 
xlabel('Número de Volums de Control'); ylabel('\alpha mitjà del tub [W/m^2·K]');
title('Evolució Coef. Transferència', 'Color', 'k');
set(ax4, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'GridColor', [0.5 0.5 0.5]);

% 
% % Gràfics de validació d'independència de malla
% figure('Name', 'Estudi de Convergència Espacial', 'NumberTitle', 'off');
% 
% subplot(2,2,1);
% plot(llista_NUM, P_sortida_malla, 'r-s', 'LineWidth', 1.5, 'MarkerFaceColor', 'r');
% grid on; xlabel('Número de Volums de Control'); ylabel('Pressió de Sortida [Pa]');
% title('Evolució Pressió');
% 
% subplot(2,2,2);
% plot(llista_NUM, H_sortida_malla, 'g-^', 'LineWidth', 1.5, 'MarkerFaceColor', 'g');
% grid on; xlabel('Número de Volums de Control'); ylabel('Entalpia de Sortida [J/kg]');
% title('Evolució Entalpia');
% 
% subplot(2,2,3);
% plot(llista_NUM, T_sortida_malla, 'b-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'b');
% grid on; xlabel('Número de Volums de Control'); ylabel('Temperatura de Sortida [ºC]');
% title('Evolució Temperatura');
% 
% subplot(2,2,4);
% plot(llista_NUM, alpha_mitja_malla, 'm-d', 'LineWidth', 1.5, 'MarkerFaceColor', 'm');
% grid on; xlabel('Número de Volums de Control'); ylabel('\alpha mitjà del tub [W/m^2·K]');
% title('Evolució Coef. Transferència');


% %% EXPORTACIÓ DELS RESULTATS A EXCEL
% 
% % 1. Creem una taula de MATLAB combinant tots els vectors de resultats
% taula_excel = table(llista_NUM', P_sortida_malla', H_sortida_malla', T_sortida_malla', alpha_mitja_malla', ...
%     'VariableNames', {'NUM_VC', 'P_Sortida_Pa', 'H_Sortida_J_kg', 'T_Sortida_C', 'alpha_mitja_W_m2K'});
% 
% % 2. Definim el nom del fitxer Excel on es guardarà
% nom_fitxer_excel = 'Estudi_Independencia_Malla.xlsx';
% 
% % 3. Exportem la taula al fitxer de forma nativa
% writetable(taula_excel, nom_fitxer_excel, 'Sheet', 'Resultats Malla');
% 
% disp(['Resultats exportats correctament a l''arxiu: ', nom_fitxer_excel]);

% %% DIBUIX DE LES CAMPANES DE SATURACIÓ (DIAGRAMES T-h I P-h)
% 
% % 1. Generació dels vectors de la campana de saturació de l'R134a
% T_crit = PropsSI('T_critical', 'P', 0, 'H', 0, Fluid); % Temperatura crítica (K)
% T_min = 200; % Límit inferior de temperatura per al gràfic (K)
% 
% T_sat_vector = linspace(T_min, T_crit - 0.1, 200); % Evitem el punt crític exacte per estabilitat
% h_l_campana = zeros(1, length(T_sat_vector));
% h_g_campana = zeros(1, length(T_sat_vector));
% P_sat_vector = zeros(1, length(T_sat_vector));
% 
% for i = 1:length(T_sat_vector)
%     P_sat_vector(i) = PropsSI('P', 'T', T_sat_vector(i), 'Q', 0, Fluid);
%     h_l_campana(i)   = PropsSI('H', 'T', T_sat_vector(i), 'Q', 0, Fluid) / 1000; % Convertim a kJ/kg
%     h_g_campana(i)   = PropsSI('H', 'T', T_sat_vector(i), 'Q', 1, Fluid) / 1000; % Convertim a kJ/kg
% end
% 
% % Unim les línies de líquid i vapor saturat per fer la forma de campana contínua
% h_campana = [h_l_campana, fliplr(h_g_campana)];
% T_campana = [T_sat_vector - 273.15, fliplr(T_sat_vector - 273.15)]; % Passat a ºC
% P_campana = [P_sat_vector / 10^5, fliplr(P_sat_vector / 10^5)];     % Passat a bar
% 
% % Convertim els vectors del teu fluid de l'última simulació a les unitats dels gràfics
% H_tub_kJ = H / 1000;      % De J/kg a kJ/kg
% P_tub_bar = P / 10^5;     % De Pa a bar
% 
% % 2. Gràfic del Diagrama Temperatura - Entalpia (T-h)
% figure('Name', 'Evolució del fluid al Diagrama T-h', 'NumberTitle', 'off');
% plot(h_campana, T_campana, 'k-', 'LineWidth', 2, 'DisplayName', 'Campana Sat. R134a');
% hold on; grid on;
% plot(H_tub_kJ, T, 'b-o', 'LineWidth', 2, 'MarkerSize', 4, 'MarkerFaceColor', 'b', 'DisplayName', 'Evolució al tub (300 VC)');
% xlabel('Entalpia específica, h [kJ/kg]');
% ylabel('Temperatura, T [ºC]');
% title('Diagrama Temperatura - Entalpia (T-h)');
% legend('Location', 'best');
% 
% % 3. Gràfic del Diagrama Pressió - Entalpia (P-h)
% figure('Name', 'Evolució del fluid al Diagrama P-h', 'NumberTitle', 'off');
% plot(h_campana, P_campana, 'k-', 'LineWidth', 2, 'DisplayName', 'Campana Sat. R134a');
% hold on; grid on;
% plot(H_tub_kJ, P_tub_bar, 'r-o', 'LineWidth', 2, 'MarkerSize', 4, 'MarkerFaceColor', 'r', 'DisplayName', 'Evolució al tub (300 VC)');
% xlabel('Entalpia específica, h [kJ/kg]');
% ylabel('Pressió, P [bar]');
% title('Diagrama Pressió - Entalpia (P-h)');
% legend('Location', 'best');