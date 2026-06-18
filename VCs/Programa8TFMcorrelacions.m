%% PROGRAMA 8 - Condicions

clear all
close all
clc

% Instal·lació de CoolProp per obtenir propietats termodinàmiques
% pyversion
% [v,e] = pyversion; system([e,' -m pip uninstall -y CoolProp']);
% [v,e] = pyversion; system([e,' -m pip install --user -U CoolProp']);

import py.CoolProp.CoolProp.*

%% DADES INICIALITZACIÓ
Fluid = 'R134a';

%% GEOMETRIA I SUPOSICIONS
% Geometria
NUM = 50;
L = 2.5;             % longitud total del tub [m]
D = 0.002;           % diàmetre interior [m]  (Kattan et al.)
S = pi * D^2/4;       % secció transversal [m²]
Perim = pi * D;       % perímetre interior [m]
dz = L / NUM;         % longitud volum de control [m]
A = Perim*dz;         % àrea mullada del VC [m²]
abs_rug = 15e-7;      % rugositat absoluta del tub [m]

% Temps molt gran (quasi estacionari)
dt = 1e9;

% Suposicions
g = 9.81;
theta = 0;            % tub horitzontal

%% VECTOR FLUX DE CALOR
Qdot = zeros(1, NUM);   % flux de calor per volum de control 

%% CONVERGÈNCIA
error_tol = 1e-5;
N = 100;



%% CONDICIÓ FRONTERA
% tipus_condicio = 1 --> flux de calor constant q [W/m²]
% tipus_condicio = 2 --> temperatura exterior fixa Text [ºC]

tipus_condicio = 1;

if tipus_condicio == 2 % Temperatura exterior fixa
    Text_ref = 210; 
    Text = Text_ref * ones(1, NUM);  % Creem el vector ple amb el valor de referència
    
    P_in = 10^6;                     % [Pa]
    %h_in = 500000;                   % J/kg
    h_in = 1100000;                  % [J/kg]
    mdot_in = 0.005;                 % [kg/s]

    % Nota: q no es defineix aquí perquè es calcularà en funció de alpha_vc
    
else
    % Flux de calor constant 
    q = 15000;                       % [W/m²]
    Text = zeros(1, NUM);            % S'omplirà a posteriori durant el bucle
    
    % Pressió de saturació a T_sat = 5 ºC
    T_sat_in = 5 + 273.15;           % [K]
    P_in = double(PropsSI('P','T',T_sat_in,'Q',0,Fluid)); 
    
    % Velocitat màssica i cabal màssic
    G = 400;                         % [kg/(m²·s)]
    mdot_in = G * S;                 % [kg/s]
    
    % Entalpia d'entrada 
    h_l_in = double(PropsSI('H','P',P_in,'Q',0,Fluid));
    h_g_in = double(PropsSI('H','P',P_in,'Q',1,Fluid));
    h_in   = h_l_in + 1000;
    
end

%% Inicialització de Vectors de Propietats
P = zeros(1, NUM+1);
H = zeros(1, NUM+1);
T = zeros(1, NUM+1);
mdot = zeros(1, NUM+1);
% Vectors de fase per a nodes i volums de control
fase    = zeros(1, NUM+1);
fase_vc = zeros(1, NUM);
alpha_vc    = zeros(1, NUM);
xg_vc_vec   = zeros(1, NUM);   % qualitat de vapor per VC (per al gràfic)

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
    cp(k) = PropsSI('C','P',P(k),'H',H(k),Fluid);
    T(k) = PropsSI('T','P',P(k),'H',H(k),Fluid) - 273.15;

    [mu(k), lambda(k)] = calc_propietats_transport(P(k), H(k), Fluid); 
    xg(k) = calc_xg(P(k), H(k), Fluid);
    fase(k) = calc_fase(P(k), H(k), Fluid);  % Càlcul de fase al node k
    
    % Fracció de buits amb correlació de Premoli (1971)
    mu_l_k  = PropsSI('V','P',P(k),'Q',0,Fluid);
    sigma_k = PropsSI('I','P',P(k),'Q',0,Fluid);
    G_k = mdot(k) / S;
    epsilon_g(k) = calc_void_premoli(G_k, D, xg(k), rho_l(k), rho_g(k), mu_l_k, sigma_k);
    rho_tp(k) = calc_rho(P(k), H(k), Fluid, G_k, D, xg(k), rho_l(k), rho_g(k), mu_l_k, sigma_k);
    v(k) = mdot(k) / (rho_tp(k) * S);
end


%% DADES INICIALITZACIÓ OLD
rho_tp_old = rho_tp;
xg_old = xg;
mdot_old = mdot;
rho_g_old = rho_g;
rho_l_old = rho_l;
epsilon_g_old = epsilon_g;
P_old = P;
H_old = H;

% Inicialitzem també els vectors variables vc_old (mida NUM=10)
rho_tp_vc_old = zeros(1, NUM);
xg_old_vc = zeros(1, NUM);
mdot_old_vc = zeros(1, NUM);
rho_g_old_vc = zeros(1, NUM);
rho_l_old_vc = zeros(1, NUM);
epsilon_g_old_vc = zeros(1, NUM);
P_old_vc = zeros(1, NUM);
h_old_vc = zeros(1, NUM);

for j = 1:NUM
    P_old_vc(j) = (P(j) + P(j+1)) / 2;
    h_old_vc(j) = (H(j) + H(j+1)) / 2;
    
    rho_l_old_vc(j) = PropsSI('D','P',P_old_vc(j),'Q',0,Fluid);
    rho_g_old_vc(j) = PropsSI('D','P',P_old_vc(j),'Q',1,Fluid);

    xg_old_vc(j) = calc_xg(P_old_vc(j), h_old_vc(j), Fluid);
    
    mdot_old_vc(j) = (mdot(j) + mdot(j+1)) / 2;

    % Fracció de buits amb correlació de Premoli (1971)
    mu_l_old_vc_j  = PropsSI('V','P',P_old_vc(j),'Q',0,Fluid);
    sigma_old_vc_j = PropsSI('I','P',P_old_vc(j),'Q',0,Fluid);
    G_vc_j = mdot_old_vc(j) / S;
    epsilon_g_old_vc(j) = calc_void_premoli(G_vc_j, D, xg_old_vc(j), rho_l_old_vc(j), rho_g_old_vc(j), mu_l_old_vc_j, sigma_old_vc_j);
    rho_tp_vc_old(j) = calc_rho(P_old_vc(j), h_old_vc(j), Fluid, G_vc_j, D, xg_old_vc(j), rho_l_old_vc(j), rho_g_old_vc(j), mu_l_old_vc_j, sigma_old_vc_j);

end




%% BUCLE ITERATIU 
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
    
        for i = 1:N    % Bucle Intern 
        iter_count = iter_count + 1;
        
        
        %% Definició de variables de convergència
        P_ite = P(j+1);
        H_ite = H(j+1);
        mdot_ite = mdot(j+1);

                %% PROPIETATS VC MEDIA 
        mdot_vc = (mdot(j) + mdot(j+1)) / 2;
        Pvc = (P(j) + P(j+1)) / 2;
        Hvc = (H(j) + H(j+1)) / 2;
        
        rho_l_vc = PropsSI('D','P',Pvc,'Q',0,Fluid);
        rho_g_vc = PropsSI('D','P',Pvc,'Q',1,Fluid);
        cp_vc = PropsSI('C','P',Pvc,'H',Hvc,Fluid);   
        [mu_vc, lambda_vc] = calc_propietats_transport(Pvc, Hvc, Fluid);
        
        Tvc = PropsSI('T','P',Pvc,'H',Hvc,Fluid) - 273.15;

        xg_vc  = calc_xg(Pvc, Hvc, Fluid);
        fase_vc(j) = calc_fase_vc(Pvc, Hvc, Fluid);  % Actualització de fase_vc al VC j
        xg_vc_vec(j) = max(0, min(1, xg_vc));        % guardem per al gràfic

        %% Propietats de líquid saturat per a bifàsic (Gungor-Winterton + Premoli)
        % Propietats en l'estat de saturació
        mu_l_vc     = PropsSI('V','P',Pvc,'Q',0,Fluid);
        mu_g_vc     = PropsSI('V','P',Pvc,'Q',1,Fluid);    % viscositat vapor saturat
        cp_l_vc     = PropsSI('C','P',Pvc,'Q',0,Fluid);
        lambda_l_vc = PropsSI('L','P',Pvc,'Q',0,Fluid);
        h_l_vc      = PropsSI('H','P',Pvc,'Q',0,Fluid);
        h_g_vc      = PropsSI('H','P',Pvc,'Q',1,Fluid);
        hfg_vc      = h_g_vc - h_l_vc;                     % salt entàlpic
        sigma_vc    = PropsSI('I','P',Pvc,'Q',0,Fluid);

        % Fracció de buits amb correlació de Premoli
        G_vc = mdot_vc / S;
        epsilon_g_vc = calc_void_premoli(G_vc, D, xg_vc, rho_l_vc, rho_g_vc, mu_l_vc, sigma_vc);
        rho_tp_vc = calc_rho(Pvc, Hvc, Fluid, G_vc, D, xg_vc, rho_l_vc, rho_g_vc, mu_l_vc, sigma_vc);
        v_vc = mdot_vc / (rho_tp_vc * S);

        % Gungor and Winterton
        alpha_vc(j) = calc_alpha_vc(Qdot(j)/A, fase_vc(j), mdot_vc, mu_vc, D, lambda_vc, cp_vc, xg_vc, rho_l_vc, rho_g_vc, mu_l_vc, cp_l_vc, lambda_l_vc, hfg_vc);
        
        % SHAH2
       % alpha_vc(j) = calc_alpha_vc_shah2(Qdot(j)/A, fase_vc(j), mdot_vc, mu_vc, D, lambda_vc, cp_vc, xg_vc, rho_l_vc, rho_g_vc, mu_l_vc, cp_l_vc, lambda_l_vc, hfg_vc);



        
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
                Text(j) = Qdot(j) / (alpha_vc(j) * A) + Tvc;
            else
                Text(j) = Tvc;
            end
        else
            % Temperatura exterior fixa (condició original)
            Qdot(j) = alpha_vc(j) * A * (Text(j) - Tvc);
        end

        % Nova entalpia al control volum de sortida
        H(j+1) = (2*Qdot(j) - mdot(j+1)*a + mdot(j)*b + ((S*dz/dt)*c)) / ((rho_tp_vc_old(j)*S*dz)/dt + mdot(j+1) + mdot(j));
        
        %% RECÀLCUL DE LA RESTA DE VARIABLES DE PROPIETATS VECTORIALS (j+1)
        rho_l(j+1) = PropsSI('D','P',P(j+1),'Q',0, Fluid);
        rho_g(j+1) = PropsSI('D','P',P(j+1),'Q',1, Fluid);     
        cp(j+1) = PropsSI('C','P',P(j+1),'H',H(j+1), Fluid);
        T(j+1) = PropsSI('T','P',P(j+1),'H',H(j+1), Fluid) - 273.15;
        [mu(j+1), lambda(j+1)] = calc_propietats_transport(P(j+1), H(j+1), Fluid);
        
        xg(j+1) = calc_xg(P(j+1), H(j+1), Fluid); 
        fase(j+1) = calc_fase(P(j+1), H(j+1), Fluid);  % Actualització de fase al node de sortida j+1

        % Fracció de buits amb correlació de Premoli 
        mu_l_jp1  = PropsSI('V','P',P(j+1),'Q',0,Fluid);
        sigma_jp1 = PropsSI('I','P',P(j+1),'Q',0,Fluid);
        G_jp1 = mdot(j+1) / S;
        epsilon_g(j+1) = calc_void_premoli(G_jp1, D, xg(j+1), rho_l(j+1), rho_g(j+1), mu_l_jp1, sigma_jp1);
        rho_tp(j+1) = calc_rho(P(j+1), H(j+1), Fluid, G_jp1, D, xg(j+1), rho_l(j+1), rho_g(j+1), mu_l_jp1, sigma_jp1);
        v(j+1) = mdot(j+1) / (rho_tp(j+1) * S);

        %% CRITERI DE CONVERGÈNCIA
        e1 = abs(mdot(j+1) - mdot_ite);
        e2 = abs(P(j+1) - P_ite);
        e3 = abs(H(j+1) - H_ite);
        
        % RELAXACIÓ
        relax = 1;
        P(j+1) = relax*P(j+1) + (1-relax)*P_ite;
        H(j+1) = relax*H(j+1) + (1-relax)*H_ite;
        mdot(j+1) = relax*mdot(j+1) + (1-relax)*mdot_ite;

        if (e1 < error_tol && e2 < error_tol && e3 < error_tol)
            break;
        end
        
    end % Fi del bucle intern (i)
    
    % % Resultados
    disp(['--- RESULTATS PER VOLUM DE CONTROL ', num2str(j), ' ---']);
    disp(['Iteracions = ', num2str(iter_count)]);
    disp(['Pin(', num2str(j), ') = ', num2str(P(j)), ' Pa']);
    disp(['Hin(', num2str(j), ') = ', num2str(H(j)), ' J/kg']);
    disp(['Tin(', num2str(j), ') = ', num2str(T(j)), ' ºC']);
    disp(['Pout(', num2str(j), ') = ', num2str(P(j+1)), ' Pa']);
    disp(['Hout(', num2str(j), ') = ', num2str(H(j+1)), ' J/kg']);
    disp(['Tout(', num2str(j), ') = ', num2str(T(j+1)), ' ºC']);
    disp(['mdot(', num2str(j), ') = ', num2str(mdot(j+1)), ' kg/s']);
    disp(['Qdot(', num2str(j), ') = ', num2str(Qdot(j)), ' W'])
    disp(['---------']);
    disp(['Pvc(', num2str(j), ') = ', num2str(Pvc), ' Pa']);
    disp(['Hvc(', num2str(j), ') = ', num2str(Hvc), ' J/kg']);
    disp(['Tvc(', num2str(j), ') = ', num2str(Tvc), ' ºC']);
    disp(['---------']);
    disp(['rho_tp_in(', num2str(j), ') = ', num2str(rho_tp(j)), ' kg/m3'])
    disp(['rho_tp_out(', num2str(j), ') = ', num2str(rho_tp(j+1)), ' kg/m3'])
    disp(['rho_g_in(', num2str(j), ') = ', num2str(rho_g(j)), ' kg/m3'])
    disp(['rho_g_out(', num2str(j), ') = ', num2str(rho_g(j+1)), ' kg/m3'])
    disp(['rho_l_in(', num2str(j), ') = ', num2str(rho_l(j)), ' kg/m3'])
    disp(['rho_l_out(', num2str(j), ') = ', num2str(rho_l(j+1)), ' kg/m3'])
    disp(['alpha_vc(', num2str(j), ') = ', num2str(alpha_vc(j)), ' '])
    disp(['tau(', num2str(j), ') = ', num2str(tau), ' ']);
    disp(['---------']);
    disp(['v_in(', num2str(j), ') = ', num2str(v(j)), ' m/s']);
    disp(['v_out(', num2str(j), ') = ', num2str(v(j+1)), ' m/s']);
    disp(['v_vc(', num2str(j), ') = ', num2str(v_vc), ' m/s']);
    disp(['---------']);
    disp(['xg_in(', num2str(j), ') = ', num2str(xg(j)), ' '])
    disp(['xg_out(', num2str(j), ') = ', num2str(xg(j+1)), ' '])
    disp(['epsilon_g_in(', num2str(j), ') = ', num2str(epsilon_g(j)), ' '])
    disp(['epsilon_g_out(', num2str(j), ') = ', num2str(epsilon_g(j+1)), ' '])
    % Resultats de fase i transferència de calor
    disp(['fase_in(', num2str(j), ')  = ', num2str(fase(j))])
    disp(['fase_out(', num2str(j), ') = ', num2str(fase(j+1))])
    disp(['fase_vc(', num2str(j), ')  = ', num2str(fase_vc(j))])
    disp(['alpha_vc(', num2str(j), ') = ', num2str(alpha_vc(j)), ' W/(m²·K)'])
    disp(['xg_vc(', num2str(j), ')   = ', num2str(xg_vc_vec(j))])
    disp(['Text(', num2str(j), ')    = ', num2str(Text(j)), ' ºC'])
    disp(['Qdot(', num2str(j), ')    = ', num2str(Qdot(j)), ' W'])
    disp(' ');

end % Fi del bucle extern (j)

%% MOSTRAR RESULTATS FINALS
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
    title(sprintf('Fluid=%s, G=%.0f kg/(m^2 s), Tsat=%.0f °C, D=%.2f mm, q=%.1f kW/m^2', Fluid, G, T_sat_in-273.15, D*1000, q/1000), 'FontSize', 12);
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

% EXPORTACIÓ A EXCEL
% Crear una taula amb les dades
ResumTaula = table((1:NUM)', xg_vc_vec', alpha_vc', 'VariableNames', {'Volum_Control', 'Qualitat_Vapor_xg', 'Alpha_W_m2K'});

% Definir el nom de l'arxiu
nom_arxiu = 'ResultatsGràfic.xlsx';

% Escriure la taula a un fitxer Excel
writetable(ResumTaula, nom_arxiu);

fprintf('\nLes dades s''han exportat correctament a l''arxiu: %s\n', nom_arxiu);