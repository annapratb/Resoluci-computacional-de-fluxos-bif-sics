%% PROGRAMA 9 - Model Transitori 
clear all
close all
clc

import py.CoolProp.CoolProp.*

%% DADES INICIALITZACIÓ
Fluid = 'R134a';

%% GEOMETRIA I SUPOSICIONS
% Geometria
NUM = 20;
L = 17.25;             % longitud total del tub [m]
D = 0.012;             % diàmetre interior [m]  (Kattan et al.)
S = pi * D^2/4;       % secció transversal [m²]
Perim = pi * D;       % perímetre interior [m]
dz = L / NUM;         % longitud volum de control [m]
A = Perim*dz;         % àrea mullada del VC [m²]
abs_rug = 15e-7;      % rugositat absoluta del tub [m]

% Suposicions
g = 9.81;
theta = 0;            % tub horitzontal

%% PARÀMETRES TEMPORALS
dt       = 250;            % Pas de temps [s]
t_max    = 1000;           % Temps màxim de simulació [s]
time     = 0;              % Comptador de temps [s]

%% CONVERGÈNCIA
error_tol = 1e-3;
N = 100;

Tsat_in = 10 + 273.15;           
P_in = double(PropsSI('P','T',Tsat_in,'Q',0,Fluid)); 
G = 300;                         
mdot_in = G * S;                 
h_l_in = double(PropsSI('H','P',P_in,'Q',0,Fluid));
h_g_in = double(PropsSI('H','P',P_in,'Q',1,Fluid));
h_in   = h_l_in + 1000;

%% CONDICIÓ FRONTERA
% tipus_condicio = 1 --> flux de calor constant q [W/m²]
% tipus_condicio = 2 --> temperatura exterior fixa Text [ºC]

tipus_condicio = 1;

if tipus_condicio == 2 % Temperatura exterior fixa
    Text_ref = 210; 
    Text = Text_ref * ones(1, NUM);  % Creem el vector ple amb el valor de referència
   
    % q no es defineix aquí perquè es calcularà en funció de alpha_vc
    
else
    % Flux de calor constant 
    q = 10000;                       % [W/m²]
    Text = zeros(1, NUM);            % S'omplirà a posteriori durant el bucle
   
end

%% Inicialització de Vectors de Propietats
P = zeros(1, NUM+1);
H = zeros(1, NUM+1);
T = zeros(1, NUM+1);
mdot = zeros(1, NUM+1);

% Vectors per a nodes i volums de control
fase    = zeros(1, NUM+1);
fase_vc = zeros(1, NUM);
alpha_vc    = zeros(1, NUM);
xg_vc_vec   = zeros(1, NUM);   % qualitat de vapor per VC (per al gràfic)
Tvc_vec = zeros(1, NUM);
Qdot = zeros(1, NUM);   % flux de calor per volum de control 

% Condicions de contorn inicials - posició 1
P(1) = P_in;
H(1) = h_in;
mdot(1) = mdot_in;

% Inicialització de la resta del vector amb les condicions d'entrada
P(2:NUM+1) = P(1);
H(2:NUM+1) = H(1);
mdot(2:NUM+1) = mdot(1);

% Inicialització de vectors de Variables de propietats (mida 11)
rho_l = zeros(1, NUM+1);
rho_g = zeros(1, NUM+1);
cp = zeros(1, NUM+1);
mu = zeros(1, NUM+1);
lambda = zeros(1, NUM+1);
epsilon_g = zeros(1, NUM+1);
xg = zeros(1, NUM+1);
rho_tp = zeros(1, NUM+1);
v = zeros(1, NUM+1);

% Càlcul de propietats d'entrada per a tot el vector
for k = 1:NUM+1
    rho_l(k) = PropsSI('D','P',P(k),'Q',0,Fluid);
    rho_g(k) = PropsSI('D','P',P(k),'Q',1,Fluid);
    h_l_k = PropsSI('H','P',P(k),'Q',0,Fluid);
    h_g_k = PropsSI('H','P',P(k),'Q',1,Fluid);
    
    T(k) = PropsSI('T','P',P(k),'H',H(k),Fluid) - 273.15;

    [mu(k), lambda(k), cp(k)] = calc_propietats_transport(P(k), H(k), Fluid);  
    xg(k) = calc_xg(P(k), H(k), Fluid, h_l_k, h_g_k);
    fase(k) = calc_fase(P(k), H(k), Fluid, h_l_k, h_g_k);  % Càlcul de fase al node k
    
    % Fracció de buits amb correlació de Premoli
    mu_l_k  = PropsSI('V','P',P(k),'Q',0,Fluid);
    sigma_k = PropsSI('I','P',P(k),'Q',0,Fluid);
    G_k = mdot(k) / S;
    epsilon_g(k) = calc_void_premoli(G_k, D, xg(k), rho_l(k), rho_g(k), mu_l_k, sigma_k);
    rho_tp(k) = calc_rho(P(k), H(k), Fluid, h_l_k, h_g_k, rho_l(k), rho_g(k), epsilon_g(k));
    v(k) = mdot(k) / (rho_tp(k) * S);
end


%% DADES INICIALITZACIÓ OLD
mdot_old = mdot;
P_old = P;
H_old = H;
rho_tp_old = rho_tp;
rho_g_old = rho_g;
rho_l_old = rho_l;
xg_old = xg;
epsilon_g_old = epsilon_g;

% Inicialitzem també els vectors variables vc_old (mida NUM=10)
mdot_old_vc = zeros(1, NUM);
P_old_vc = zeros(1, NUM);
h_old_vc = zeros(1, NUM);
rho_tp_vc_old = zeros(1, NUM);
rho_g_old_vc = zeros(1, NUM);
rho_l_old_vc = zeros(1, NUM);
xg_old_vc = zeros(1, NUM);
epsilon_g_old_vc = zeros(1, NUM);

for j = 1:NUM
    mdot_old_vc(j) = (mdot(j) + mdot(j+1)) / 2;
    P_old_vc(j) = (P(j) + P(j+1)) / 2;
    h_old_vc(j) = (H(j) + H(j+1)) / 2;
    
    rho_l_old_vc(j) = PropsSI('D','P',P_old_vc(j),'Q',0,Fluid);
    rho_g_old_vc(j) = PropsSI('D','P',P_old_vc(j),'Q',1,Fluid);
    h_l_old_vc_j = PropsSI('H','P',P_old_vc(j),'Q',0,Fluid);
    h_g_old_vc_j = PropsSI('H','P',P_old_vc(j),'Q',1,Fluid);
    xg_old_vc(j) = calc_xg(P_old_vc(j), h_old_vc(j), Fluid, h_l_old_vc_j, h_g_old_vc_j);
    % Fracció de buits amb correlació de Premoli 
    mu_l_old_vc_j  = PropsSI('V','P',P_old_vc(j),'Q',0,Fluid);
    sigma_old_vc_j = PropsSI('I','P',P_old_vc(j),'Q',0,Fluid);
    G_vc_j = mdot_old_vc(j) / S;
    epsilon_g_old_vc(j) = calc_void_premoli(G_vc_j, D, xg_old_vc(j), rho_l_old_vc(j), rho_g_old_vc(j), mu_l_old_vc_j, sigma_old_vc_j);
    rho_tp_vc_old(j) = calc_rho(P_old_vc(j), h_old_vc(j), Fluid, h_l_old_vc_j, h_g_old_vc_j, rho_l_old_vc(j), rho_g_old_vc(j), epsilon_g_old_vc(j));

end


%% Inicialització historial temporal (pel gràfic)
max_hist = min(round(t_max/dt), 50000); % màxim passos a desar. Calcula quants punts de dades tindrem i ho talla a 50.000
hist_time = nan(1, max_hist);  % Crea una fila de "buits" (NaN significa Not a Number).
hist_mdot_out = nan(1, max_hist);  % mdot al node de sortida
hist_P_out = nan(1, max_hist);  % P al node de sortida
hist_H_out = nan(1, max_hist);  % H al node de sortida
hist_err = nan(1, max_hist);  % error màxim global per pas
n_hist = 0; % comptador de passos desats. De moment, encara no s'ha escrit res a la "llibreta".

%% BUCLE TEMPORAL

for t_step = 1:(t_max/dt)   % Bucle temporal fins a t_max/dt passos
    % Avançar en el temps
    time = time + dt;

    % Guardar estat de l'instant anterior
    mdot_prev = mdot_old;
    P_prev = P_old;
    H_prev = H_old;
    
    for j = 1:NUM      % Bucle Extern
    iter_count = 0;

    % Valors per a la primera iteració
    P(j+1)=P(j);
    H(j+1)=H(j);
    xg(j+1)=xg(j);
    rho_l(j+1)=rho_l(j);
    rho_g(j+1)=rho_g(j);
    epsilon_g(j+1)=epsilon_g(j);
    alpha_vc(j)=1000;
    
    %% PROPIETATS VC MEDIA 
    mdot_vc = (mdot(j) + mdot(j+1)) / 2;
    Pvc = (P(j) + P(j+1)) / 2;
    Hvc = (H(j) + H(j+1)) / 2;
    
    rho_l_vc = PropsSI('D','P',Pvc,'Q',0,Fluid);
    rho_g_vc = PropsSI('D','P',Pvc,'Q',1,Fluid);
    h_l_vc = PropsSI('H','P',Pvc,'Q',0,Fluid);
    h_g_vc = PropsSI('H','P',Pvc,'Q',1,Fluid);
    hfg_vc = h_g_vc - h_l_vc;                     % salt entàlpic

    [mu_vc, lambda_vc, cp_vc] = calc_propietats_transport(Pvc, Hvc, Fluid);

    Tvc = PropsSI('T','P',Pvc,'H',Hvc,Fluid) - 273.15;
    Tvc_vec(j) = Tvc;  % Guardem la temperatura del VC per a l'acoblament amb la paret

    xg_vc  = calc_xg(Pvc, Hvc, Fluid, h_l_vc, h_g_vc);
    fase_vc(j) = calc_fase_vc(Pvc, Hvc, Fluid, h_l_vc, h_g_vc);  % Actualització de fase_vc al VC j
    xg_vc_vec(j) = max(0, min(1, xg_vc));        % guardem per al gràfic

    % Propietats de líquid saturat per a bifàsic (Gungor-Winterton + Premoli)
    % Propietats en l'estat de saturació
    mu_l_vc = PropsSI('V','P',Pvc,'Q',0,Fluid);
    mu_g_vc = PropsSI('V','P',Pvc,'Q',1,Fluid);    % viscositat vapor saturat
    cp_l_vc = PropsSI('C','P',Pvc,'Q',0,Fluid);
    lambda_l_vc = PropsSI('L','P',Pvc,'Q',0,Fluid);
    sigma_vc = PropsSI('I','P',Pvc,'Q',0,Fluid);

    % Fracció de buits amb correlació de Premoli
    G_vc = mdot_vc / S;
    epsilon_g_vc = calc_void_premoli(G_vc, D, xg_vc, rho_l_vc, rho_g_vc, mu_l_vc, sigma_vc);
    rho_tp_vc = calc_rho(Pvc, Hvc, Fluid, h_l_vc, h_g_vc, rho_l_vc, rho_g_vc, epsilon_g_vc);
    v_vc = mdot_vc / (rho_tp_vc * S);
    
    %% Càlcul alpha
    % Defineix el mètode: 'Gungor' o 'SHAH'
    correlacio_calcul = 'GungorWinterton'; 

    if strcmp(correlacio_calcul, 'GungorWinterton')
    % Gungor and Winterton
    alpha_vc(j) = calc_alpha_vc(Qdot(j)/A, fase_vc(j), mdot_vc, mu_vc, D, lambda_vc, cp_vc, xg_vc, rho_l_vc, rho_g_vc, mu_l_vc, cp_l_vc, lambda_l_vc, hfg_vc);

    elseif strcmp(correlacio_calcul, 'SHAH')
    % SHAH
    alpha_vc(j) = calc_alpha_vc_shah2(Qdot(j)/A, fase_vc(j), mdot_vc, mu_vc, D, lambda_vc, cp_vc, xg_vc, rho_l_vc, rho_g_vc, mu_l_vc, cp_l_vc, lambda_l_vc, hfg_vc);

    end

    for i = 1:N    % Bucle Intern 
        iter_count = iter_count + 1;
        
        %% Definició de variables de convergència
        P_ite = P(j+1);
        H_ite = H(j+1);
        mdot_ite = mdot(j+1);
        
        %% EQUACIÓ DE LA MASSA
        mdot(j+1) = mdot(j) - ((rho_tp(j) - rho_tp_old(j))/dt) * S * dz;
        
        %% EQUACIÓ MOMENTUM
        tau = calc_tau(fase_vc(j), fase(j), mdot_vc, D, S, abs_rug, rho_tp_vc, mu_vc, xg_vc, rho_l_vc, rho_g_vc, mu_l_vc, mu_g_vc, sigma_vc);
        
        % TERME 1: convectiu sortida
        term_1 = (mdot(j+1)/S)^2 * ( xg(j+1)^2/(rho_g(j+1)*epsilon_g(j+1)) + (1-xg(j+1))^2/(rho_l(j+1)*(1-epsilon_g(j+1))) );

        % TERME 2: convectiu entrada
        term_2 = (mdot(j)/S)^2 * ( xg(j)^2/(rho_g(j)*epsilon_g(j)) + (1-xg(j))^2/(rho_l(j)*(1-epsilon_g(j))) );
        
        % TERME 3: transitori
        term_t = (mdot_vc - mdot_old_vc(j))/(S*dt);

        term_3 = term_t + (tau * Perim / S) + rho_tp_vc * g * sin(theta);
        
        % Equació final
        P(j+1) = P(j) - term_1 + term_2 - term_3 * dz;

        %% EQUACIÓ ENERGIA
        % TERME a: 
        a = ( xg(j+1) / 2 ) * ( (mdot(j+1)^2 * xg(j+1)^2) / (rho_g(j+1)^2 * epsilon_g(j+1)^2 * S^2) - (mdot(j)^2 * xg(j)^2) / (rho_g(j)^2 * epsilon_g(j)^2 * S^2) ) + ...
            (( 1 - xg(j+1)) / 2) * ( (mdot(j+1)^2 * (1 - xg(j+1))^2) / (rho_l(j+1)^2 * (1 - epsilon_g(j+1))^2 * S^2) - (mdot(j)^2 * (1 - xg(j))^2) / (rho_l(j)^2 * ( 1 - epsilon_g(j) )^2 * S^2) ) ...
            + (g*sin(theta)*dz - H(j)); 
            
        % TERME b: 
        b = ( xg(j) / 2) * ( (mdot(j)^2 * xg(j)^2) / (rho_g(j)^2 * epsilon_g(j)^2 * S^2) - (mdot(j+1)^2 * xg(j+1)^2) / (rho_g(j+1)^2 * epsilon_g(j+1)^2 * S^2) ) + ...
            ((1 - xg(j)) / 2) * ( (mdot(j)^2 * (1 - xg(j))^2) / (rho_l(j)^2 * (1 - epsilon_g(j))^2 * S^2) - (mdot(j+1)^2 * (1 - xg(j+1))^2) / (rho_l(j+1)^2 * (1 - epsilon_g(j+1) )^2 * S^2) ) ...
            + (-g*sin(theta)*dz + H(j));

        % TERME c: transitori 
        c = rho_tp_vc_old(j) * (H_old(j+1) + H_old(j) - H(j)) + ( (P(j) + P(j+1)) - (P_old(j) + P_old(j+1)) ) ...
            - (rho_g_old_vc(j) * epsilon_g_old_vc(j))/2 * ( (mdot(j)^2  * xg(j)^2) / (rho_g(j)^2 * epsilon_g(j)^2 * S^2) + (mdot(j+1)^2 * xg(j+1)^2) / (rho_g(j+1)^2 * epsilon_g(j+1)^2 * S^2) ) ...
            - (rho_l_old_vc(j) * (1 - epsilon_g_old_vc(j)))/2 * ( (mdot(j)^2  * (1 - xg(j))^2) / (rho_l(j)^2 * (1 - epsilon_g(j))^2 * S^2) + (mdot(j+1)^2 * (1 - xg(j+1))^2) / (rho_l(j+1)^2 * (1 - epsilon_g(j+1))^2 * S^2) ) ...
            + (rho_g_old_vc(j) * epsilon_g_old_vc(j))/2 * ( (mdot_old(j)^2 * xg_old(j)^2) / (rho_g_old(j)^2 * epsilon_g_old(j)^2 * S^2) + (mdot_old(j+1)^2 * xg_old(j+1)^2) / (rho_g_old(j+1)^2 * epsilon_g_old(j+1)^2 * S^2) ) ...
            + (rho_l_old_vc(j) * (1 - epsilon_g_old_vc(j)))/2 * ( (mdot_old(j)^2 * (1 - xg_old(j))^2) / (rho_l_old(j)^2 * (1 - epsilon_g_old(j))^2 * S^2) + (mdot_old(j+1)^2 * (1 - xg_old(j+1))^2) / (rho_l_old(j+1)^2 * (1 - epsilon_g_old(j+1))^2 * S^2) );

        %% FLUX DE CALOR i CONDICIÓ DE FRONTERA
        if tipus_condicio == 1
            % Flux de calor constant
            Qdot(j) = q * A;
            % Temperatura exterior equivalent
            if alpha_vc(j) > 0
                Text(j) = Qdot(j) / (alpha_vc(j) * A) + Tvc_vec(j);
            else
                Text(j) = Tvc_vec(j);
            end
        else
            % Temperatura exterior fixa (condició original)
            Qdot(j) = alpha_vc(j) * A * (Text(j) - Tvc_vec(j));
        end

        % Nova entalpia al control volum de sortida
        H(j+1) = (2*Qdot(j) - mdot(j+1)*a + mdot(j)*b + ((S*dz/dt)*c)) / ((rho_tp_vc_old(j)*S*dz)/dt + mdot(j+1) + mdot(j));
        

        %% RECÀLCUL DE LA RESTA DE VARIABLES DE PROPIETATS VECTORIALS (j+1)
        rho_l(j+1) = PropsSI('D','P',P(j+1),'Q',0, Fluid);
        rho_g(j+1) = PropsSI('D','P',P(j+1),'Q',1, Fluid);     
        h_l_jp1 = PropsSI('H','P',P(j+1),'Q',0, Fluid);
        h_g_jp1 = PropsSI('H','P',P(j+1),'Q',1, Fluid);
        
        [mu(j+1), lambda(j+1), cp(j+1)] = calc_propietats_transport(P(j+1), H(j+1), Fluid);
                
        xg(j+1) = calc_xg(P(j+1), H(j+1), Fluid, h_l_jp1, h_g_jp1); 

        % Fracció de buits amb correlació de Premoli 
        mu_l_jp1  = PropsSI('V','P',P(j+1),'Q',0,Fluid);
        sigma_jp1 = PropsSI('I','P',P(j+1),'Q',0,Fluid);
        G_jp1 = mdot(j+1) / S;
        epsilon_g(j+1) = calc_void_premoli(G_jp1, D, xg(j+1), rho_l(j+1), rho_g(j+1), mu_l_jp1, sigma_jp1);
        rho_tp(j+1) = calc_rho(P(j+1), H(j+1), Fluid, h_l_jp1, h_g_jp1, rho_l(j+1), rho_g(j+1), epsilon_g(j+1));
        v(j+1) = mdot(j+1) / (rho_tp(j+1) * S);

        %% CRITERI DE CONVERGÈNCIA
        e1 = abs(mdot(j+1) - mdot_ite);
        e2 = abs(P(j+1) - P_ite);
        e3 = abs(H(j+1) - H_ite);
        
        % RELAXACIÓ
        relax = 0.5;
        P(j+1) = relax*P(j+1) + (1-relax)*P_ite;
        H(j+1) = relax*H(j+1) + (1-relax)*H_ite;
        mdot(j+1) = relax*mdot(j+1) + (1-relax)*mdot_ite;

        if (e1 < error_tol && e2 < error_tol && e3 < error_tol)
            break;
        end
        
    end % Fi del bucle intern (i)

    %% PROPIETATS POST-CONVERGÈNCIA (no necessàries en cada iteració)
    cp(j+1) = PropsSI('C','P',P(j+1),'H',H(j+1), Fluid);
    T(j+1) = PropsSI('T','P',P(j+1),'H',H(j+1), Fluid) - 273.15;
    h_l_jp1 = PropsSI('H','P',P(j+1),'Q',0, Fluid);
    h_g_jp1 = PropsSI('H','P',P(j+1),'Q',1, Fluid);
    fase(j+1) = calc_fase(P(j+1), H(j+1), Fluid, h_l_jp1, h_g_jp1); 
    
    % Resultats
    disp(['--- RESULTATS PER VOLUM DE CONTROL ', num2str(j), ' ---']);
    disp(['Iteracions = ', num2str(iter_count)]);
    disp(['mdot_in(', num2str(j), ') = ', num2str(mdot(j)), ' kg/s']);
    disp(['Pin(', num2str(j), ') = ', num2str(P(j)), ' Pa']);
    disp(['Hin(', num2str(j), ') = ', num2str(H(j)), ' J/kg']);
    disp(['Tin(', num2str(j), ') = ', num2str(T(j)), ' ºC']);
    disp(['Pout(', num2str(j), ') = ', num2str(P(j+1)), ' Pa']);
    disp(['Hout(', num2str(j), ') = ', num2str(H(j+1)), ' J/kg']);
    disp(['Tout(', num2str(j), ') = ', num2str(T(j+1)), ' ºC']);
    disp(['mdot_out(', num2str(j), ') = ', num2str(mdot(j+1)), ' kg/s']);
    % disp(['Qdot(', num2str(j), ') = ', num2str(Qdot(j)), ' W'])
    % disp(['---------']);
    % disp(['Pvc(', num2str(j), ') = ', num2str(Pvc), ' Pa']);
    % disp(['Hvc(', num2str(j), ') = ', num2str(Hvc), ' J/kg']);
    % disp(['Tvc(', num2str(j), ') = ', num2str(Tvc), ' ºC']);
    % disp(['---------']);
    % disp(['rho_tp_in(', num2str(j), ') = ', num2str(rho_tp(j)), ' kg/m3'])
    % disp(['rho_tp_out(', num2str(j), ') = ', num2str(rho_tp(j+1)), ' kg/m3'])
    % disp(['rho_g_in(', num2str(j), ') = ', num2str(rho_g(j)), ' kg/m3'])
    % disp(['rho_g_out(', num2str(j), ') = ', num2str(rho_g(j+1)), ' kg/m3'])
    % disp(['rho_l_in(', num2str(j), ') = ', num2str(rho_l(j)), ' kg/m3'])
    % disp(['rho_l_out(', num2str(j), ') = ', num2str(rho_l(j+1)), ' kg/m3'])
    % disp(['alpha_vc(', num2str(j), ') = ', num2str(alpha_vc(j)), ' '])
    % disp(['tau(', num2str(j), ') = ', num2str(tau), ' ']);
    % disp(['---------']);
    % disp(['v_in(', num2str(j), ') = ', num2str(v(j)), ' m/s']);
    % disp(['v_out(', num2str(j), ') = ', num2str(v(j+1)), ' m/s']);
    % disp(['v_vc(', num2str(j), ') = ', num2str(v_vc), ' m/s']);
    % disp(['---------']);
    disp(['xg_in(', num2str(j), ') = ', num2str(xg(j)), ' '])
    disp(['xg_out(', num2str(j), ') = ', num2str(xg(j+1)), ' '])
    disp(['xg_vc(', num2str(j), ')   = ', num2str(xg_vc_vec(j))])
    disp(['epsilon_g_in(', num2str(j), ') = ', num2str(epsilon_g(j)), ' '])
    disp(['epsilon_g_out(', num2str(j), ') = ', num2str(epsilon_g(j+1)), ' '])
    disp(['fase_in(', num2str(j), ')  = ', num2str(fase(j))])
    disp(['fase_out(', num2str(j), ') = ', num2str(fase(j+1))])
    disp(['fase_vc(', num2str(j), ')  = ', num2str(fase_vc(j))])
    disp(['alpha_vc(', num2str(j), ') = ', num2str(alpha_vc(j)), ' W/(m²·K)'])
    disp(['Text(', num2str(j), ')    = ', num2str(Text(j)), ' ºC'])
    disp(['Qdot(', num2str(j), ')    = ', num2str(Qdot(j)), ' W'])
    disp(' ');

    end % Fi del bucle espacial (j)

    %% CRITERI DE CONVERGÈNCIA TEMPORAL
    err_mdot = max(abs(mdot - mdot_prev));
    err_P = max(abs(P - P_prev));
    err_H = max(abs(H - H_prev));

    %% Desar l'estat actual a l'historial (pel gràfic) 
    % Omple els buits amb els resultats reals de cada segon
    if n_hist < max_hist
        n_hist = n_hist + 1;
        hist_time(n_hist) = time;
        hist_mdot_out(n_hist) = mdot(end); % l'"end" agafa automàticament l'últim node del tub
        hist_P_out(n_hist) = P(end);
        hist_H_out(n_hist) = H(end);
        hist_err(n_hist) = max([err_mdot, err_P, err_H]); % busca quin dels tres errors és el més gran en aquell instant i el guarda.
    end

    % Informe periòdic cada 1 pas de temps: mod(..., 1): calcula el residu d'una divisió. Dividir per 1 i mirar si el residu és 0 és com: "Escriu només quan el número sigui rodó".
    % Si el divisor és 1, escriurà a cada pas de temps. Si fos 10, escriuria cada 10 passos de temps
    if mod(round(time/dt), 1) == 0         % time/dt: Calcula en quin número de pas de temps està
        fprintf('t = %8.1f s | err_mdot = %.3e | err_P = %.3e | err_H = %.3e\n', ...
                time, err_mdot, err_P, err_H);
    end

    % Comprovació de l'equilibri estacionari
    if (err_mdot <= error_tol && err_P <= error_tol && err_H <= error_tol)
        fprintf('\n>>> Estat estacionari assolit en t = %.2f s (%.0f passos de temps)\n', time, time/dt);
        break;
    end

    %% ACTUALITZACIÓ DE VARIABLES OLD
    % Variables als nodes (mida NUM+1)
    mdot_old = mdot;
    P_old = P;
    H_old = H;
    rho_tp_old = rho_tp;
    rho_g_old = rho_g;
    rho_l_old = rho_l;
    xg_old = xg;
    epsilon_g_old = epsilon_g;

    % Recàlcul de propietats als centres dels old_vc
    for j = 1:NUM
        mdot_old_vc(j)  = (mdot(j) + mdot(j+1)) / 2;
        P_old_vc(j) = (P(j) + P(j+1)) / 2;
        h_old_vc(j) = (H(j) + H(j+1)) / 2;

        rho_l_old_vc(j) = PropsSI('D','P',P_old_vc(j),'Q',0,Fluid);
        rho_g_old_vc(j) = PropsSI('D','P',P_old_vc(j),'Q',1,Fluid);
        h_l_old_vc_j = PropsSI('H','P',P_old_vc(j),'Q',0,Fluid);
        h_g_old_vc_j = PropsSI('H','P',P_old_vc(j),'Q',1,Fluid);
        xg_old_vc(j) = calc_xg(P_old_vc(j), h_old_vc(j), Fluid, h_l_old_vc_j, h_g_old_vc_j);        
        mu_l_old_vc_j  = PropsSI('V','P',P_old_vc(j),'Q',0,Fluid);
        sigma_old_vc_j = PropsSI('I','P',P_old_vc(j),'Q',0,Fluid);
        G_vc_j = mdot_old_vc(j) / S;
        epsilon_g_old_vc(j) = calc_void_premoli(G_vc_j, D, xg_old_vc(j), rho_l_old_vc(j), rho_g_old_vc(j), mu_l_old_vc_j, sigma_old_vc_j);
        rho_tp_vc_old(j) = calc_rho(P_old_vc(j), h_old_vc(j), Fluid, h_l_old_vc_j, h_g_old_vc_j, rho_l_old_vc(j), rho_g_old_vc(j), epsilon_g_old_vc(j));
    end

end % Fi del bucle temporal 


%% MOSTRAR RESULTATS FINALS (t = temps d'equilibri)
fprintf('\n=== RESULTATS FINALS (t = %.2f s) ===\n', time);
disp('============= RESULTATS FINALS =============');
for k = 1:NUM+1
    fprintf('Node %2d: P = %6.0f Pa, H = %8.0f J/kg, T = %5.2f ºC, fase = %2.0f, mdot = %5.4f kg/s\n', k, P(k), H(k), T(k), fase(k), mdot(k));
end

%% GRÀFIC: Coeficient de Transferència de Calor vs Qualitat de Vapor
figure('Name', 'Validació', 'NumberTitle', 'off');
plot(xg_vc_vec, alpha_vc, 'b-o', 'LineWidth', 2, 'MarkerSize', 7, 'MarkerFaceColor', 'b');
xlabel('Vapor quality [-]', 'FontSize', 13);
ylabel('Heat transfer coefficient [W/(m^2 K)]', 'FontSize', 13);
if tipus_condicio == 2 
   title(sprintf('Fluid = %s, D = %.2f mm', Fluid, D), 'FontSize', 12);  
else
    title(sprintf('Fluid=%s, G=%.0f kg/(m^2 s), Tsat=%.0f °C, D=%.2f mm, q=%.1f kW/m^2', Fluid, G, Tsat_in-273.15, D*1000, q/1000), 'FontSize', 12);
end
% Límits eixos
grid on;
if tipus_condicio == 2 
    set(gca, 'FontSize', 12, 'XLim', [0 1], 'YLim', [0 125000]);
else
    set(gca, 'FontSize', 12, 'XLim', [0 1], 'YLim', [0 8000]); % Eix X de 0 a 1, Eix Y de 0 a 8000
end

fprintf('\nResum per VC: xg_vc i alpha_vc\n');
for j = 1:NUM
    fprintf('  VC %2d: xg = %.4f,  alpha = %.2f W/(m²·K)\n', j, xg_vc_vec(j), alpha_vc(j));
end

%% EXPORTACIÓ A EXCEL
% % Crear una taula amb les dades
% ResumTaula = table((1:NUM)', xg_vc_vec', alpha_vc', 'VariableNames', {'Volum_Control', 'Qualitat_Vapor_xg', 'Alpha_W_m2K'});
% 
% % Definir el nom de l'arxiu
% nom_arxiu = 'ResultatsGràfic.xlsx';
% 
% % Escriure la taula a un fitxer Excel
% writetable(ResumTaula, nom_arxiu);
% 
% fprintf('\nLes dades s''han exportat correctament a l''arxiu: %s\n', nom_arxiu);


%% GRÀFIC: Evolució temporal de mdot, P i H fins a l'equilibri estacionari

if n_hist > 0
    t_plot   = hist_time(1:n_hist);
    m_plot   = hist_mdot_out(1:n_hist);
    P_plot   = hist_P_out(1:n_hist);
    H_plot   = hist_H_out(1:n_hist);
    err_plot = hist_err(1:n_hist);

    figure('Name', 'Evolució Temporal fins a l''Equilibri', 'NumberTitle', 'off', ...
           'Position', [50 50 950 750]);

    % Subgràfic 1: Cabal màssic de sortida 
    subplot(4,1,1);       % 4 files, 1 columna, gràfic 1
    plot(t_plot, m_plot * 1000, 'b-', 'LineWidth', 1.5);
    ylabel('mdot_{out} [g/s]', 'FontSize', 11);
    title(sprintf('Evolució temporal — Fluid: %s | dt = %g s | Equilibri a t = %.1f s', ...
                  Fluid, dt, time), 'FontSize', 11);
    grid on; box on;
    set(gca, 'FontSize', 10);

    % Subgràfic 2: Pressió de sortida
    subplot(4,1,2);      % 4 files, 1 columna, gràfic 2
    plot(t_plot, P_plot / 1e5, 'r-', 'LineWidth', 1.5);
    ylabel('P_{out} [bar]', 'FontSize', 11);
    grid on; box on;
    set(gca, 'FontSize', 10);

    % Subgràfic 3: Entalpia de sortida
    subplot(4,1,3);      % 4 files, 1 columna, gràfic 3
    plot(t_plot, H_plot / 1e3, 'Color', [0.1 0.6 0.1], 'LineWidth', 1.5); % color verd fosc personalitzat
    ylabel('H_{out} [kJ/kg]', 'FontSize', 11);
    grid on; box on;
    set(gca, 'FontSize', 10);

    % Subgràfic 4: Error de convergència (escala logarítmica)
    subplot(4,1,4);     % 4 files, 1 columna, gràfic 3
    semilogy(t_plot, err_plot, 'y-', 'LineWidth', 1);  % eix vertical en escala logaritmica (l'error). Línia negra continua.
    hold on;
    semilogy([t_plot(1) t_plot(end)], [error_tol error_tol], 'r--', 'LineWidth', 1.5); % discontinua vermella
    ylabel('Error max. [-]', 'FontSize', 11);
    xlabel('Temps [s]', 'FontSize', 11);
    legend('Error global', sprintf('Tolerancia = %.0e', error_tol), ...
           'Location', 'northeast', 'FontSize', 9);
    grid on; box on;
    set(gca, 'FontSize', 10);

    % Marcar el punt d'equilibri en tots els subgràfics
    t_eq = time;
    for isp = 1:4     % Entra a cada subgràfic 
        subplot(4,1,isp);   % Activa cada subgràfic  per afegir el punt d'equilibri
        xline(t_eq, 'm--', sprintf('t_{eq} = %.1f s', t_eq), 'LabelVerticalAlignment', 'bottom', 'FontSize', 9, 'LineWidth', 1);  % Dibuixa una línia vertical (x line) en el valor de temps t_eq
    end

end

