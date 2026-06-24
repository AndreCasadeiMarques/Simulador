%% =====================================================================
% SIMULADOR GNC PARA eVTOL (GD-350)
% Descrição: Laço principal de simulação GNC (Guiamento, Navegação e 
% Controle) - Padrão Arquitetural Estrito (IMAV-M) - Zero Cálculos na Main
% =====================================================================

%% CONFIGURAÇÃO DO AMBIENTE E PARÂMETROS
% Limpa o workspace e carrega os parâmetros físicos e aerodinâmicos do veículo
clear; clc; close all;

% Adiciona todas as subpastas ao Path do MATLAB
addpath(genpath(pwd));

% Gera e carrega as structs de parâmetros por domínio
run('parameters.m');
load('parameters.mat', 'sMav', 'sControl', 'sGuidance', 'sSim');


%% INSTANCIAÇÃO DE OBJETOS (INJEÇÃO DE DEPENDÊNCIA)
% Inicializa os módulos de software injetando apenas o seu domínio de dados
mav     = cMav(sMav);
ctrl    = cControl(sControl);
guide   = cGuidance(sGuidance);
plotter = cPlotter();


%% CONDIÇÃO INICIAL (TRIMAGEM DE HOVER - DELEGADA À CLASSE DE PLANTA)
eta_hover = mav.getHoverInput();

% Inicializa rotores 1-8 em hover, e rotores 9-10 (horizontais) em zero
mav.varpi = [eta_hover * sMav.km(1) * ones(8,1); 0; 0];
eta_prev  = [eta_hover * ones(8,1); 0; 0];
delta_aero = zeros(3,1);                         % Deflexão de superfícies de controle [rad]


%% CONFIGURAÇÃO DO TEMPO E PRÉ-ALOCAÇÃO DE MEMÓRIA (LOGGING PASSIVO)
t = 0:sSim.Ts:sSim.t_sim;
N = length(t);

% --- Estados da Planta e Referências ---
hist.r       = zeros(3, N);     % Posição real [m]
hist.r_bar   = zeros(3, N);     % Posição de referência (Guiamento) [m]
hist.v       = zeros(3, N);     % Velocidade linear [m/s]
hist.v_bar   = zeros(3, N);     % Velocidade de referência [m/s]
hist.a_bar   = zeros(3, N);     % Aceleração de referência [m/s^2]
hist.q       = zeros(4, N);     % Atitude em quatérnio (escalar no final)
hist.q_bar   = zeros(4, N);     % Atitude de referência comandada (escalar no final)
hist.w       = zeros(3, N);     % Velocidade angular [rad/s]

% --- Atuadores e Controle ---
hist.varpi   = zeros(sSim.n_r, N); % Rotação real dos motores [rad/s]
hist.eta     = zeros(sSim.n_r, N); % Comando enviado aos motores [0 a 1]
hist.f_cmd   = zeros(3, N);     % Força resultante comandada [N]
hist.tau_cmd = zeros(3, N);     % Torque resultante comandado [N.m]
hist.f_star  = zeros(sSim.n_r, N); % Empuxo ideal exigido por rotor [N]

% --- Diagnóstico Aerodinâmico e Físico ---
hist.alpha   = zeros(1, N);     % Ângulo de ataque [rad]
hist.beta    = zeros(1, N);     % Ângulo de derrapagem [rad]
hist.pdin    = zeros(1, N);     % Pressão dinâmica [Pa]
hist.f_real  = zeros(3, N);     % Força fisicamente realizada [N]
hist.tau_real= zeros(3, N);     % Torque fisicamente realizado [N.m]


%% LAÇO PRINCIPAL DE SIMULAÇÃO (DISCRETE-TIME LOOP - MAESTRO)
disp('>> Simulando...');

for k = 1:N
    % 1. Log dos Estados Atuais (Feedback e Logging Passivos)
    hist.r(:, k)     = mav.r;
    hist.v(:, k)     = mav.v;
    hist.q(:, k)     = mav.q; 
    hist.w(:, k)     = mav.w;
    hist.varpi(:, k) = mav.varpi;
    
    % 2. Atualização de Perturbações e Esforços Externos (Aerodinâmica)
    mav.updateDisturbances(eta_prev, delta_aero);
    
    % 3. Guiamento (Geração de Trajetória)
    ref = guide.getCommand(mav.r, mav.v);
    hist.r_bar(:, k) = ref.r_bar;
    hist.v_bar(:, k) = ref.v_bar;
    hist.a_bar(:, k) = ref.a_bar;
    
    % 4. Controle (Cálculo de Ação e Alocação) - PLANTA DESACOPLADA DO CONTROLADOR
    [eta, log_ctrl] = ctrl.compute(ref, mav.r, mav.v, mav.q, mav.w, mav.f_aero, mav.tau_aero);
    
    % 5. Log de Ações do Controlador (Logging Passivo)
    hist.eta(:, k)     = eta; 
    hist.q_bar(:, k)   = D2q(log_ctrl.D_cmd); % Conversão DCM -> Quatérnio
    hist.f_cmd(:, k)   = log_ctrl.f_cmd;
    hist.tau_cmd(:, k) = log_ctrl.tau_cmd;
    hist.f_star(:, k)  = log_ctrl.f_star;
    
    % 6. Planta (Integração RK4 - Propagação da Física)
    mav.integrate(eta, delta_aero);
    eta_prev = eta; % Armazena o comando para causalidade do motor no próximo passo
    
    % 7. Diagnóstico Aerodinâmico e Auditoria Física (Logging Passivo de Caixas Pretas)
    diag_data = mav.getDiagnostics();
    hist.alpha(1, k)    = diag_data.alpha;
    hist.beta(1, k)     = diag_data.beta;
    hist.pdin(1, k)     = diag_data.pdin;
    hist.f_real(:, k)   = diag_data.f_real;
    hist.tau_real(:, k) = diag_data.tau_real;
end


%% PÓS-PROCESSAMENTO E EXPORTAÇÃO
% Salva os resultados para análise e aciona o módulo de plotagem
save('Outputs/sim_data.mat', 'hist', 't', 'sMav', 'sControl', 'sGuidance', 'sSim');
disp('>> Simulação concluída. Dados salvos em Outputs/sim_data.mat');

% Compila matriz para exportação externa (ex: Python, Excel)
dados_exportacao = [t', hist.r', hist.r_bar', hist.v', hist.q', hist.q_bar', ...
                    hist.w', hist.f_cmd', hist.tau_cmd', hist.f_star', ...
                    hist.varpi', hist.eta'];
                
writematrix(dados_exportacao, 'Outputs/sim_data.csv');

% Geração de Gráficos (Delegação Total de Plotagem)
plotter.plotAll(t, hist);