% =========================================================================
% SCRIPT DE MATLAB: Interpolació amb eliminació de duplicats (xg = 1)
% =========================================================================

% 1. Configurar els noms dels arxius
arxiu_original = 'ResultatsGràfic.xlsx'; % El teu arxiu nou
arxiu_interp = 'dades_interpolades_v2.xlsx'; % L'arxiu que es crearà

% 2. Llegir les dades de l'Excel
taula_original = readtable(arxiu_original);

% 3. NETEJA DE DUPLICATS: Ens quedem només amb els punts on xg és únic
% 'stable' fa que es mantingui l'ordre original de les files.
% Això agafarà el primer xg = 1 (fila 34) i ignorarà els de sota (35 a 50).
[~, indexs_unics] = unique(taula_original.Qualitat_Vapor_xg, 'stable');

xg_neteja = taula_original.Qualitat_Vapor_xg(indexs_unics);
alpha_neteja = taula_original.Alpha_W_m2K(indexs_unics);

% 4. Crear el nou vector xg d'interès (de 0 en 0.05 fins a 1)
xg_nou = (0:0.05:1)';

% 5. Realitzar la interpolació lineal amb les dades netes
% S'utilitza 'extrap' per calcular el valor inicial a xg = 0
alpha_nou = interp1(xg_neteja, alpha_neteja, xg_nou, 'linear', 'extrap');

% 6. Crear la nova taula amb els resultats nets i ben estructurats
taula_nova = table(xg_nou, alpha_nou, ...
    'VariableNames', {'Qualitat_Vapor_xg', 'Alpha_W_m2K'});

% 7. Guardar els resultats en un nou arxiu d'Excel
writetable(taula_nova, arxiu_interp);

% Missatge de confirmació a la consola
fprintf('Procés finalitzat correctament!\n');
fprintf('S''han eliminat els valors duplicats de xg=1 per poder interpolar.\n');
fprintf('Arxiu desat com a: %s\n', arxiu_interp);