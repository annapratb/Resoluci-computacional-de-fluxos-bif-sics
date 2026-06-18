%% PROGRAMA 9 - Model Transitori (Estudi d'independència de dt)
clear all
close all
clc
import py.CoolProp.CoolProp.*

%% DADES INICIALITZACIÓ
Fluid = 'R134a';

%% GEOMETRIA I SUPOSICIONS
NUM = 50;
L = 17.25;             % longitud total del tub [m]
D = 0.012;             % diàmetre interior [m]
S = pi * D^2/4;        % secció transversal [m²]
Perim = pi * D;        % perímetre interior [m]
dz = L / NUM;          % longitud volum de control [m]
A = Perim*dz;          % àrea mullada del VC [m²]
abs_rug = 15e-7;       % rugositat absoluta del tub [m]
g = 9.81;
theta = 0;             % tub horitzontal

%% PARÀMETRES TEMPORALS DE L'ESTUDI
% Definim la llista de passos de temps (dt) que volem comparar
% Com més gran sigui el dt, més es desviarà el perfil transitori de la realitat.
llista_dt = [500, 250, 100, 50, 10, 5, 3];  
%llista_dt = [500, 250, 100, 50];  
t_max    = 700;           % Temps màxim de simulació [s]

%% CONVERGÈNCIA
error_tol = 1e-5;
N = 100;
Tsat_in = 10 + 273.15;           
P_in = double(PropsSI('P','T',Tsat_in,'Q',0,Fluid)); 
G = 300;                         
mdot_in = G * S;                 
h_l_in = double(PropsSI('H','P',P_in,'Q',0,Fluid));
h_g_in = double(PropsSI('H','P',P_in,'Q',1,Fluid));
h_in   = h_l_in + 1000;

%% CONDICIÓ FRONTERA
tipus_condicio = 1;
if tipus_condicio == 2
    Text_ref = 210; 
    Text = Text_ref * ones(1, NUM);  
else
    q = 10000;                       % [W/m²]
    Text = zeros(1, NUM);            
end

%% STRUCT PER GUARDAR TOTS ELS RESULTATS DELS DIFERENTS DT
% Aquí guardarem l'historial de cada simulació per poder-les comparar al final
resultats_dt = struct();

%% BUCLE PRINCIPAL PER CADA PAS DE TEMPS (dt)
for idx_dt = 1:length(llista_dt)
    dt = llista_dt(idx_dt);
    fprintf('\n=======================================\n');
    fprintf(' INICIANT SIMULACIÓ AMB dt = %g s\n', dt);
    fprintf('=======================================\n');
    
    time = 0;
    
    % Inicialització de Vectors de Propietats per a cada simu
    P = zeros(1, NUM+1); H = zeros(1, NUM+1); T = zeros(1, NUM+1); mdot = zeros(1, NUM+1);
    fase = zeros(1, NUM+1); fase_vc = zeros(1, NUM); alpha_vc = zeros(1, NUM);
    xg_vc_vec = zeros(1, NUM); Tvc_vec = zeros(1, NUM); Qdot = zeros(1, NUM);   
    
    P(1) = P_in; H(1) = h_in; mdot(1) = mdot_in;
    P(2:NUM+1) = P(1); H(2:NUM+1) = H(1); mdot(2:NUM+1) = mdot(1);
    
    rho_l = zeros(1, NUM+1); rho_g = zeros(1, NUM+1); cp = zeros(1, NUM+1);
    mu = zeros(1, NUM+1); lambda = zeros(1, NUM+1); epsilon_g = zeros(1, NUM+1);
    xg = zeros(1, NUM+1); rho_tp = zeros(1, NUM+1); v = zeros(1, NUM+1);
    
    for k = 1:NUM+1
        rho_l(k) = PropsSI('D','P',P(k),'Q',0,Fluid);
        rho_g(k) = PropsSI('D','P',P(k),'Q',1,Fluid);
        h_l_k = PropsSI('H','P',P(k),'Q',0,Fluid); h_g_k = PropsSI('H','P',P(k),'Q',1,Fluid);
        T(k) = PropsSI('T','P',P(k),'H',H(k),Fluid) - 273.15;
        [mu(k), lambda(k), cp(k)] = calc_propietats_transport(P(k), H(k), Fluid);  
        xg(k) = calc_xg(P(k), H(k), Fluid, h_l_k, h_g_k);
        fase(k) = calc_fase(P(k), H(k), Fluid, h_l_k, h_g_k);  
        mu_l_k  = PropsSI('V','P',P(k),'Q',0,Fluid); sigma_k = PropsSI('I','P',P(k),'Q',0,Fluid);
        G_k = mdot(k) / S;
        epsilon_g(k) = calc_void_premoli(G_k, D, xg(k), rho_l(k), rho_g(k), mu_l_k, sigma_k);
        rho_tp(k) = calc_rho(P(k), H(k), Fluid, h_l_k, h_g_k, rho_l(k), rho_g(k), epsilon_g(k));
        v(k) = mdot(k) / (rho_tp(k) * S);
    end
    
    % Estat OLD inicial
    mdot_old = mdot; P_old = P; H_old = H; rho_tp_old = rho_tp; rho_g_old = rho_g; rho_l_old = rho_l; xg_old = xg; epsilon_g_old = epsilon_g;
    mdot_old_vc = zeros(1, NUM); P_old_vc = zeros(1, NUM); h_old_vc = zeros(1, NUM); rho_tp_vc_old = zeros(1, NUM); rho_g_old_vc = zeros(1, NUM); rho_l_old_vc = zeros(1, NUM); xg_old_vc = zeros(1, NUM); epsilon_g_old_vc = zeros(1, NUM);
    for j = 1:NUM
        mdot_old_vc(j) = (mdot(j) + mdot(j+1)) / 2;
        P_old_vc(j) = (P(j) + P(j+1)) / 2;
        h_old_vc(j) = (H(j) + H(j+1)) / 2;
        rho_l_old_vc(j) = PropsSI('D','P',P_old_vc(j),'Q',0,Fluid);
        rho_g_old_vc(j) = PropsSI('D','P',P_old_vc(j),'Q',1,Fluid);
        h_l_old_vc_j = PropsSI('H','P',P_old_vc(j),'Q',0,Fluid); h_g_old_vc_j = PropsSI('H','P',P_old_vc(j),'Q',1,Fluid);
        xg_old_vc(j) = calc_xg(P_old_vc(j), h_old_vc(j), Fluid, h_l_old_vc_j, h_g_old_vc_j);
        mu_l_old_vc_j  = PropsSI('V','P',P_old_vc(j),'Q',0,Fluid); sigma_old_vc_j = PropsSI('I','P',P_old_vc(j),'Q',0,Fluid);
        G_vc_j = mdot_old_vc(j) / S;
        epsilon_g_old_vc(j) = calc_void_premoli(G_vc_j, D, xg_old_vc(j), rho_l_old_vc(j), rho_g_old_vc(j), mu_l_old_vc_j, sigma_old_vc_j);
        rho_tp_vc_old(j) = calc_rho(P_old_vc(j), h_old_vc(j), Fluid, h_l_old_vc_j, h_g_old_vc_j, rho_l_old_vc(j), rho_g_old_vc(j), epsilon_g_old_vc(j));
    end
    
    % Inicialització historial temporal per aquest dt específic
    max_hist = round(t_max/dt) + 1;
    hist_time = nan(1, max_hist);  
    hist_P_out = nan(1, max_hist);  
    hist_H_out = nan(1, max_hist);  
    hist_T_out = nan(1, max_hist);
    hist_alpha_mean = nan(1, max_hist); % Guardarem la mitjana del tub d'alpha per simplificar el gràfic temporal
    n_hist = 0;
    
    %% BUCLE TEMPORAL
    for t_step = 1:(t_max/dt)   
        time = time + dt;
        mdot_prev = mdot_old; P_prev = P_old; H_prev = H_old;
        
        for j = 1:NUM      % Bucle Espacial
            iter_count = 0;
            P(j+1)=P(j); H(j+1)=H(j); xg(j+1)=xg(j); rho_l(j+1)=rho_l(j); rho_g(j+1)=rho_g(j); epsilon_g(j+1)=epsilon_g(j); alpha_vc(j)=1000;
            
            mdot_vc = (mdot(j) + mdot(j+1)) / 2; Pvc = (P(j) + P(j+1)) / 2; Hvc = (H(j) + H(j+1)) / 2;
            rho_l_vc = PropsSI('D','P',Pvc,'Q',0,Fluid); rho_g_vc = PropsSI('D','P',Pvc,'Q',1,Fluid);
            h_l_vc = PropsSI('H','P',Pvc,'Q',0,Fluid); h_g_vc = PropsSI('H','P',Pvc,'Q',1,Fluid); hfg_vc = h_g_vc - h_l_vc;                     
            [mu_vc, lambda_vc, cp_vc] = calc_propietats_transport(Pvc, Hvc, Fluid);
            Tvc = PropsSI('T','P',Pvc,'H',Hvc,Fluid) - 273.15; Tvc_vec(j) = Tvc;  
            xg_vc  = calc_xg(Pvc, Hvc, Fluid, h_l_vc, h_g_vc);
            fase_vc(j) = calc_fase_vc(Pvc, Hvc, Fluid, h_l_vc, h_g_vc);  
            xg_vc_vec(j) = max(0, min(1, xg_vc));        
            
            mu_l_vc = PropsSI('V','P',Pvc,'Q',0,Fluid); mu_g_vc = PropsSI('V','P',Pvc,'Q',1,Fluid);    
            cp_l_vc = PropsSI('C','P',Pvc,'Q',0,Fluid); lambda_l_vc = PropsSI('L','P',Pvc,'Q',0,Fluid); sigma_vc = PropsSI('I','P',Pvc,'Q',0,Fluid);
            G_vc = mdot_vc / S;
            epsilon_g_vc = calc_void_premoli(G_vc, D, xg_vc, rho_l_vc, rho_g_vc, mu_l_vc, sigma_vc);
            rho_tp_vc = calc_rho(Pvc, Hvc, Fluid, h_l_vc, h_g_vc, rho_l_vc, rho_g_vc, epsilon_g_vc);
            
            correlacio_calcul = 'GungorWinterton'; 
            if strcmp(correlacio_calcul, 'GungorWinterton')
                alpha_vc(j) = calc_alpha_vc(Qdot(j)/A, fase_vc(j), mdot_vc, mu_vc, D, lambda_vc, cp_vc, xg_vc, rho_l_vc, rho_g_vc, mu_l_vc, cp_l_vc, lambda_l_vc, hfg_vc);
            elseif strcmp(correlacio_calcul, 'SHAH')
                alpha_vc(j) = calc_alpha_vc_shah2(Qdot(j)/A, fase_vc(j), mdot_vc, mu_vc, D, lambda_vc, cp_vc, xg_vc, rho_l_vc, rho_g_vc, mu_l_vc, cp_l_vc, lambda_l_vc, hfg_vc);
            end
            
            for i = 1:N    % Bucle Intern Convergència
                iter_count = iter_count + 1;
                P_ite = P(j+1); H_ite = H(j+1); mdot_ite = mdot(j+1);
                
                mdot(j+1) = mdot(j) - ((rho_tp(j) - rho_tp_old(j))/dt) * S * dz;
                tau = calc_tau(fase_vc(j), fase(j), mdot_vc, D, S, abs_rug, rho_tp_vc, mu_vc, xg_vc, rho_l_vc, rho_g_vc, mu_l_vc, mu_g_vc, sigma_vc);
                term_1 = (mdot(j+1)/S)^2 * ( xg(j+1)^2/(rho_g(j+1)*epsilon_g(j+1)) + (1-xg(j+1))^2/(rho_l(j+1)*(1-epsilon_g(j+1))) );
                term_2 = (mdot(j)/S)^2 * ( xg(j)^2/(rho_g(j)*epsilon_g(j)) + (1-xg(j))^2/(rho_l(j)*(1-epsilon_g(j))) );
                term_t = (mdot_vc - mdot_old_vc(j))/(S*dt);
                term_3 = term_t + (tau * Perim / S) + rho_tp_vc * g * sin(theta);
                P(j+1) = P(j) - term_1 + term_2 - term_3 * dz;
                
                a = ( xg(j+1) / 2 ) * ( (mdot(j+1)^2 * xg(j+1)^2) / (rho_g(j+1)^2 * epsilon_g(j+1)^2 * S^2) - (mdot(j)^2 * xg(j)^2) / (rho_g(j)^2 * epsilon_g(j)^2 * S^2) ) + ...
                    (( 1 - xg(j+1)) / 2) * ( (mdot(j+1)^2 * (1 - xg(j+1))^2) / (rho_l(j+1)^2 * (1 - epsilon_g(j+1))^2 * S^2) - (mdot(j)^2 * (1 - xg(j))^2) / (rho_l(j)^2 * ( 1 - epsilon_g(j) )^2 * S^2) ) + (g*sin(theta)*dz - H(j)); 
                b = ( xg(j) / 2) * ( (mdot(j)^2 * xg(j)^2) / (rho_g(j)^2 * epsilon_g(j)^2 * S^2) - (mdot(j+1)^2 * xg(j+1)^2) / (rho_g(j+1)^2 * epsilon_g(j+1)^2 * S^2) ) + ...
                    ((1 - xg(j)) / 2) * ( (mdot(j)^2 * (1 - xg(j))^2) / (rho_l(j)^2 * (1 - epsilon_g(j))^2 * S^2) - (mdot(j+1)^2 * xg(j+1)^2) / (rho_g(j+1)^2 * (1 - epsilon_g(j+1) )^2 * S^2) ) + (-g*sin(theta)*dz + H(j));
                c = rho_tp_vc_old(j) * (H_old(j+1) + H_old(j) - H(j)) + ( (P(j) + P(j+1)) - (P_old(j) + P_old(j+1)) ) ...
                    - (rho_g_old_vc(j) * epsilon_g_old_vc(j))/2 * ( (mdot(j)^2  * xg(j)^2) / (rho_g(j)^2 * epsilon_g(j)^2 * S^2) + (mdot(j+1)^2 * xg(j+1)^2) / (rho_g(j+1)^2 * epsilon_g(j+1)^2 * S^2) ) ...
                    - (rho_l_old_vc(j) * (1 - epsilon_g_old_vc(j)))/2 * ( (mdot(j)^2  * (1 - xg(j))^2) / (rho_l(j)^2 * (1 - epsilon_g(j))^2 * S^2) + (mdot(j+1)^2 * (1 - xg(j+1))^2) / (rho_l(j+1)^2 * (1 - epsilon_g(j+1))^2 * S^2) ) ...
                    + (rho_g_old_vc(j) * epsilon_g_old_vc(j))/2 * ( (mdot_old(j)^2 * xg_old(j)^2) / (rho_g_old(j)^2 * epsilon_g_old(j)^2 * S^2) + (mdot_old(j+1)^2 * xg_old(j+1)^2) / (rho_g_old(j+1)^2 * epsilon_g_old(j+1)^2 * S^2) ) ...
                    + (rho_l_old_vc(j) * (1 - epsilon_g_old_vc(j)))/2 * ( (mdot_old(j)^2 * (1 - xg_old(j))^2) / (rho_l_old(j)^2 * (1 - epsilon_g_old(j))^2 * S^2) + (mdot_old(j+1)^2 * (1 - xg_old(j+1))^2) / (rho_l_old(j+1)^2 * (1 - epsilon_g_old(j+1))^2 * S^2) );
                
                if tipus_condicio == 1
                    Qdot(j) = q * A;
                else
                    Qdot(j) = alpha_vc(j) * A * (Text(j) - Tvc_vec(j));
                end
                H(j+1) = (2*Qdot(j) - mdot(j+1)*a + mdot(j)*b + ((S*dz/dt)*c)) / ((rho_tp_vc_old(j)*S*dz)/dt + mdot(j+1) + mdot(j));
                
                rho_l(j+1) = PropsSI('D','P',P(j+1),'Q',0, Fluid); rho_g(j+1) = PropsSI('D','P',P(j+1),'Q',1, Fluid);     
                h_l_jp1 = PropsSI('H','P',P(j+1),'Q',0, Fluid); h_g_jp1 = PropsSI('H','P',P(j+1),'Q',1, Fluid);
                [mu(j+1), lambda(j+1), cp(j+1)] = calc_propietats_transport(P(j+1), H(j+1), Fluid);
                xg(j+1) = calc_xg(P(j+1), H(j+1), Fluid, h_l_jp1, h_g_jp1); 
                mu_l_jp1  = PropsSI('V','P',P(j+1),'Q',0,Fluid); sigma_jp1 = PropsSI('I','P',P(j+1),'Q',0,Fluid);
                G_jp1 = mdot(j+1) / S;
                epsilon_g(j+1) = calc_void_premoli(G_jp1, D, xg(j+1), rho_l(j+1), rho_g(j+1), mu_l_jp1, sigma_jp1);
                rho_tp(j+1) = calc_rho(P(j+1), H(j+1), Fluid, h_l_jp1, h_g_jp1, rho_l(j+1), rho_g(j+1), epsilon_g(j+1));
                
                e1 = abs(mdot(j+1) - mdot_ite); e2 = abs(P(j+1) - P_ite); e3 = abs(H(j+1) - H_ite);
                relax = 0.5;
                P(j+1) = relax*P(j+1) + (1-relax)*P_ite; H(j+1) = relax*H(j+1) + (1-relax)*H_ite; mdot(j+1) = relax*mdot(j+1) + (1-relax)*mdot_ite;
                if (e1 < error_tol && e2 < error_tol && e3 < error_tol), break; end
            end 
            T(j+1) = PropsSI('T','P',P(j+1),'H',H(j+1), Fluid) - 273.15;
            h_l_jp1 = PropsSI('H','P',P(j+1),'Q',0, Fluid); h_g_jp1 = PropsSI('H','P',P(j+1),'Q',1, Fluid);
            fase(j+1) = calc_fase(P(j+1), H(j+1), Fluid, h_l_jp1, h_g_jp1); 
        end % Fi del bucle espacial (j)
        
        %% CRITERI DE CONVERGÈNCIA TEMPORAL
        err_mdot = max(abs(mdot - mdot_prev)); err_P = max(abs(P - P_prev)); err_H = max(abs(H - H_prev));
        
        %% DESAR ESTAT ACTUAL A L'HISTORIAL
        n_hist = n_hist + 1;
        hist_time(n_hist) = time;
        hist_P_out(n_hist) = P(end);         % Pressió a la sortida del tub
        hist_H_out(n_hist) = H(end);         % Entalpia a la sortida del tub
        hist_T_out(n_hist) = T(end);         % Temperatura a la sortida del tub
        hist_alpha_mean(n_hist) = mean(alpha_vc); % Coeficient alfa mitjà del tub
        
        if (err_mdot <= error_tol && err_P <= error_tol && err_H <= error_tol)
            fprintf('>>> Estat estacionari assolit per a dt = %g s en t = %.1f s\n', dt, time);
            break;
        end
        
        % Actualització d'estats transitoris OLD
        mdot_old = mdot; P_old = P; H_old = H; rho_tp_old = rho_tp; rho_g_old = rho_g; rho_l_old = rho_l; xg_old = xg; epsilon_g_old = epsilon_g;
        for j = 1:NUM
            mdot_old_vc(j)  = (mdot(j) + mdot(j+1)) / 2;
            P_old_vc(j) = (P(j) + P(j+1)) / 2;
            h_old_vc(j) = (H(j) + H(j+1)) / 2;
            rho_l_old_vc(j) = PropsSI('D','P',P_old_vc(j),'Q',0,Fluid); rho_g_old_vc(j) = PropsSI('D','P',P_old_vc(j),'Q',1,Fluid);
            h_l_old_vc_j = PropsSI('H','P',P_old_vc(j),'Q',0,Fluid); h_g_old_vc_j = PropsSI('H','P',P_old_vc(j),'Q',1,Fluid);
            xg_old_vc(j) = calc_xg(P_old_vc(j), h_old_vc(j), Fluid, h_l_old_vc_j, h_g_old_vc_j);        
            mu_l_old_vc_j  = PropsSI('V','P',P_old_vc(j),'Q',0,Fluid); sigma_old_vc_j = PropsSI('I','P',P_old_vc(j),'Q',0,Fluid);
            G_vc_j = mdot_old_vc(j) / S;
            epsilon_g_old_vc(j) = calc_void_premoli(G_vc_j, D, xg_old_vc(j), rho_l_old_vc(j), rho_g_old_vc(j), mu_l_old_vc_j, sigma_old_vc_j);
            rho_tp_vc_old(j) = calc_rho(P_old_vc(j), h_old_vc(j), Fluid, h_l_old_vc_j, h_g_old_vc_j, rho_l_old_vc(j), rho_g_old_vc(j), epsilon_g_old_vc(j));
        end
    end % Fi del bucle temporal
    
    % Netegem els NaNs si ha convergit abans d'hora
    hist_time = hist_time(1:n_hist);
    hist_P_out = hist_P_out(1:n_hist);
    hist_H_out = hist_H_out(1:n_hist);
    hist_T_out = hist_T_out(1:n_hist);
    hist_alpha_mean = hist_alpha_mean(1:n_hist);
    
    % Guardem les dades d'aquest dt a la macro-estructura
    camp = sprintf('dt_%g', dt);
    resultats_dt.(camp).time = hist_time;
    resultats_dt.(camp).P = hist_P_out;
    resultats_dt.(camp).H = hist_H_out;
    resultats_dt.(camp).T = hist_T_out;
    resultats_dt.(camp).alpha = hist_alpha_mean;
    
end % Fi del bucle de llista_dt


%% GENERACIÓ DELS GRÀFICS COMPARATIUS (Independència del Pas de Temps - FONS BLANC)
% 1. Creem la finestra forçant el fons general a blanc ('Color', 'w')
fig = figure('Name', 'Estudi d''Independència del Pas de Temps (dt)', 'NumberTitle', 'off', ...
    'Position', [100 50 1000 800], 'Color', 'w');

colors = {'r', 'b', 'g', 'm', 'c'}; % Colors per diferenciar els diferents dt

% --- Subgràfic 1: Pressió de Sortida vs Temps ---
ax1 = subplot(4,1,1); hold on; box on; grid on;
for idx_dt = 1:length(llista_dt)
    camp = sprintf('dt_%g', llista_dt(idx_dt));
    plot(resultats_dt.(camp).time, resultats_dt.(camp).P / 1e5, colors{mod(idx_dt-1,5)+1}, 'LineWidth', 1.5);
end
ylabel('P_{out} [bar]', 'Color', 'k');
title('Evolució temporal de les variables principals segons la mida de dt', 'Color', 'k');
% Ajustos de color per a l'eix 1
set(ax1, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'GridColor', [0.5 0.5 0.5]);

% --- Subgràfic 2: Entalpia de Sortida vs Temps ---
ax2 = subplot(4,1,2); hold on; box on; grid on;
for idx_dt = 1:length(llista_dt)
    camp = sprintf('dt_%g', llista_dt(idx_dt));
    plot(resultats_dt.(camp).time, resultats_dt.(camp).H / 1e3, colors{mod(idx_dt-1,5)+1}, 'LineWidth', 1.5);
end
ylabel('H_{out} [kJ/kg]', 'Color', 'k');
% Ajustos de color per a l'eix 2
set(ax2, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'GridColor', [0.5 0.5 0.5]);

% --- Subgràfic 3: Temperatura de Sortida vs Temps ---
ax3 = subplot(4,1,3); hold on; box on; grid on;
for idx_dt = 1:length(llista_dt)
    camp = sprintf('dt_%g', llista_dt(idx_dt));
    plot(resultats_dt.(camp).time, resultats_dt.(camp).T, colors{mod(idx_dt-1,5)+1}, 'LineWidth', 1.5);
end
ylabel('T_{out} [ºC]', 'Color', 'k');
% Ajustos de color per a l'eix 3
set(ax3, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'GridColor', [0.5 0.5 0.5]);

% --- Subgràfic 4: Alpha mitjà vs Temps ---
ax4 = subplot(4,1,4); hold on; box on; grid on;
etiquetes_llegenda = cell(1, length(llista_dt));
for idx_dt = 1:length(llista_dt)
    camp = sprintf('dt_%g', llista_dt(idx_dt));
    plot(resultats_dt.(camp).time, resultats_dt.(camp).alpha, colors{mod(idx_dt-1,5)+1}, 'LineWidth', 1.5);
    etiquetes_llegenda{idx_dt} = sprintf('dt = %g s', llista_dt(idx_dt));
end
ylabel('\alpha_{mitjà} [W/m^2K]', 'Color', 'k');
xlabel('Temps [s]', 'Color', 'k');
% Ajustos de color per a l'eix 4
set(ax4, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'GridColor', [0.5 0.5 0.5]);

% Afegim la llegenda forçant que el seu fons també sigui blanc i text negre
lgd = legend(etiquetes_llegenda, 'Location', 'best');
set(lgd, 'Color', 'w', 'TextColor', 'k', 'EdgeColor', [0.3 0.3 0.3]);

% 
% 
% %% GENERACIÓ DELS GRÀFICS COMPARATIUS (Independència del Pas de Temps)
% figure('Name', 'Estudi d''Independència del Pas de Temps (dt)', 'NumberTitle', 'off', 'Position', [100 50 1000 800]);
% colors = {'r', 'b', 'g', 'm', 'c'}; % Colors per diferenciar els diferents dt
% 
% % Subgràfic 1: Pressió de Sortida vs Temps
% subplot(4,1,1); hold on; box on; grid on;
% for idx_dt = 1:length(llista_dt)
%     camp = sprintf('dt_%g', llista_dt(idx_dt));
%     plot(resultats_dt.(camp).time, resultats_dt.(camp).P / 1e5, colors{mod(idx_dt-1,5)+1}, 'LineWidth', 1.5);
% end
% ylabel('P_{out} [bar]');
% title('Evolució temporal de les variables principals segons la mida de dt');
% 
% % Subgràfic 2: Entalpia de Sortida vs Temps
% subplot(4,1,2); hold on; box on; grid on;
% for idx_dt = 1:length(llista_dt)
%     camp = sprintf('dt_%g', llista_dt(idx_dt));
%     plot(resultats_dt.(camp).time, resultats_dt.(camp).H / 1e3, colors{mod(idx_dt-1,5)+1}, 'LineWidth', 1.5);
% end
% ylabel('H_{out} [kJ/kg]');
% 
% % Subgràfic 3: Temperatura de Sortida vs Temps
% subplot(4,1,3); hold on; box on; grid on;
% for idx_dt = 1:length(llista_dt)
%     camp = sprintf('dt_%g', llista_dt(idx_dt));
%     plot(resultats_dt.(camp).time, resultats_dt.(camp).T, colors{mod(idx_dt-1,5)+1}, 'LineWidth', 1.5);
% end
% ylabel('T_{out} [ºC]');
% 
% % Subgràfic 4: Alpha mitjà vs Temps
% subplot(4,1,4); hold on; box on; grid on;
% etiquetes_llegenda = cell(1, length(llista_dt));
% for idx_dt = 1:length(llista_dt)
%     camp = sprintf('dt_%g', llista_dt(idx_dt));
%     plot(resultats_dt.(camp).time, resultats_dt.(camp).alpha, colors{mod(idx_dt-1,5)+1}, 'LineWidth', 1.5);
%     etiquetes_llegenda{idx_dt} = sprintf('dt = %g s', llista_dt(idx_dt));
% end
% ylabel('\alpha_{mitjà} [W/m^2K]');
% xlabel('Temps [s]');
% 
% % Afegim la llegenda només a l'últim subgràfic per no carregar el dibuix
% legend(etiquetes_llegenda, 'Location', 'best');
% 
% %% EXPORTACIÓ DE RESULTATS A EXCEL
% nom_fitxer = 'Estudi_Independencia_dt.xlsx';
% 
% % Si el fitxer ja existeix, l'esborrem per no barrejar dades velles
% if exist(nom_fitxer, 'file') == 2
%     delete(nom_fitxer); 
% end
% 
% fprintf('\n=======================================\n');
% fprintf(' EXPORTANT DADES A EXCEL...\n');
% fprintf('=======================================\n');
% 
% for idx_dt = 1:length(llista_dt)
%     dt_actual = llista_dt(idx_dt);
%     camp = sprintf('dt_%g', dt_actual);
% 
%     % 1. Extraiem les dades de l'estructura i les posem en format columna
%     temps_col  = resultats_dt.(camp).time';
%     P_out_col  = (resultats_dt.(camp).P / 1e5)';      % Convertit a bar
%     H_out_col  = (resultats_dt.(camp).H / 1e3)';      % Convertit a kJ/kg
%     T_out_col  = resultats_dt.(camp).T';
%     alpha_col  = resultats_dt.(camp).alpha';
% 
%     % 2. Creem una taula de MATLAB per a aquest dt específic
%     T_excel = table(temps_col, P_out_col, H_out_col, T_out_col, alpha_col, ...
%         'VariableNames', {'Temps_s', 'P_out_bar', 'H_out_kJ_kg', 'T_out_C', 'Alpha_mitja_W_m2K'});
% 
%     % 3. Guardem la taula en una pestanya de l'Excel amb el nom del dt
%     nom_pestanya = sprintf('dt = %g s', dt_actual);
%     writetable(T_excel, nom_fitxer, 'Sheet', nom_pestanya);
% 
%     fprintf('Pestanya "%s" guardada correctament.\n', nom_pestanya);
% end
% 
% fprintf('Done! El fitxer "%s" s''ha creat a la teva carpeta de treball.\n', nom_fitxer);
% fprintf('=======================================\n');