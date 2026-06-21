% =========================================================================
% TESTE UNITÁRIO DA CLASSE DE CONTROLE
% =========================================================================
clc; clear; close all;

disp('>> Iniciando Testes Unitários do Controlador...');

% 1. Carrega os parâmetros (garanta que seu parameters.m roda aqui)
parameters; 

% 2. Instancia o controlador
ctrl = cControl(p);

% 3. Cria structs vazias e estáticas para simular a aeronave e a referência
mav = struct('r', [0;0;0], 'v', [0;0;0], 'q', [1;0;0;0], 'w', [0;0;0], ...
             'f_aero', [0;0;0], 'tau_aero', [0;0;0], 'tau_tilde', [0;0;0]);
ref = struct('r_bar', [0;0;0], 'v_bar', [0;0;0], 'a_bar', [0;0;0], ...
             'alpha_bar', [0;0;0], 'w_bar', [0;0;0], 'w_dot_bar', [0;0;0]);

% =========================================================================
% TESTE 1: VOO PAIRADO (HOVER)
% Drone na origem. Referência na origem.
% =========================================================================
disp(' '); disp('--- TESTE 1: HOVER (EQUILÍBRIO) ---');
[eta, log_ctrl] = ctrl.compute(ref, mav);

fprintf('Empuxo Vertical (Tz_req): %.2f N\n', log_ctrl.f_cmd(3));
disp('Torques Comandados [tx; ty; tz]:'); disp(log_ctrl.tau_cmd');
disp('>> AVALIAÇÃO: O torque DEVE ser [0 0 0]. O empuxo DEVE ser exatamente o Peso do veículo (m*g). Se o empuxo der zero, sua compensação de gravidade sumiu.');

% =========================================================================
% TESTE 2: DEGRAU DE ALTITUDE
% Drone no Z = 0. Referência pede para subir para Z = -10 (ou +10).
% =========================================================================
disp(' '); disp('--- TESTE 2: DEGRAU DE ALTITUDE ---');
ref.r_bar = [0; 0; 10]; % Ajuste o sinal do 10 dependendo do seu referencial (Z p/ baixo = negativo)
[eta, log_ctrl] = ctrl.compute(ref, mav);

fprintf('Empuxo Vertical (Tz_req): %.2f N\n', log_ctrl.f_cmd(3));
disp('>> AVALIAÇÃO: O empuxo DEVE ser MAIOR (em módulo) que o Peso (m*g), pois o drone precisa subir. Se for menor, o sinal do erro de posição Z no código está invertido!');

% =========================================================================
% TESTE 3: ERRO DE POSIÇÃO FRONTAL (Eixo X)
% Drone na origem. Referência no X = 10.
% =========================================================================
disp(' '); disp('--- TESTE 3: DESLOCAMENTO FRONTAL ---');
ref.r_bar = [10; 0; 0]; % Quer ir para frente
[~, log_ctrl] = ctrl.compute(ref, mav);

% Converter a matriz de atitude comandada para ângulos de Euler para conferir
D_cmd = log_ctrl.D_cmd;
pitch_cmd = -asin(max(min(D_cmd(1,3), 1), -1)) * (180/pi);

fprintf('Comando de Pitch (Arfagem): %.2f graus\n', pitch_cmd);
disp('>> AVALIAÇÃO: Para ir para a frente (X+), o drone DEVE inclinar o nariz para BAIXO. O Pitch comandado não pode ser zero, deve ter um valor sensato (ex: -5 a -20 graus). Se der NaN, sua malha de atitude explodiu na extração do erro.');

% =========================================================================
% TESTE 4: ERRO PURAMENTE DE ATITUDE (PERTURBAÇÃO)
% Drone sofreu um vento e rolou 30 graus (Erro de Roll). 
% =========================================================================
disp(' '); disp('--- TESTE 4: RECUPERAÇÃO DE ROLAGEM ---');
ref.r_bar = [0;0;0]; % Reseta referência para origem
% Simulando que o drone está inclinado 30 graus em Roll (X) usando quatérnio
ang = 30 * (pi/180);
mav.q = [cos(ang/2); sin(ang/2); 0; 0]; 
[~, log_ctrl] = ctrl.compute(ref, mav);

disp('Torques Comandados [tx; ty; tz]:'); disp(log_ctrl.tau_cmd');
disp('>> AVALIAÇÃO: O torque em X (tx) DEVE ser gigantesco e com SINAL OPOSTO à perturbação para forçar o drone a voltar a zero graus. Se tx = 0 ou tx acompanhar o mesmo sinal da perturbação, seu controlador está injetando energia ao invés de amortecer (Feedback Positivo).');
disp('Empuxo comandado em cada rotor (f_star):');
disp(log_ctrl.f_star');