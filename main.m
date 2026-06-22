%% =====================================================================
% SIMULADOR GNC PARA eVTOL (GD-350)
% Descrição: Laço principal de simulação GNC (Guiamento, Navegação e 
% Controle) utilizando integração de Runge-Kutta de 4ª Ordem.
% =====================================================================

%% CONFIGURAÇÃO DO AMBIENTE E PARÂMETROS
% Limpa o workspace e carrega os parâmetros físicos e aerodinâmicos do veículo
clear; clc; close all;

% Adiciona todas as subpastas ao Path do MATLAB
addpath(genpath(pwd));

% Gera e carrega a struct de parâmetros (p)
run('parameters.m');
load('parameters.mat', 'p');


%% INSTANCIAÇÃO DE OBJETOS (INJEÇÃO DE DEPENDÊNCIA)
% Inicializa os módulos de software da arquitetura do simulador
mav     = cMav(p);
ctrl    = cControl(p);
guide   = cGuidance(p);
plotter = cPlotter();


%% CONDIÇÃO INICIAL (TRIMAGEM DE HOVER)
% Calcula analiticamente o esforço para sustentar o veículo no ar e 
% inicializa as velocidades angulares dos rotores verticais.
f_hover   = (p.m * p.g) / (8 * cos(3 * pi/180)); % Empuxo por rotor [N]
w_hover   = sqrt(f_hover / p.kf(1));             % Rotação de hover [rad/s]
eta_hover = w_hover / p.km(1);                   % Comando do motor de hover [0 a 1]

% Inicializa rotores 1-8 em hover, e rotores 9-10 (horizontais) em zero
mav.varpi = [w_hover * ones(8,1); 0; 0];
eta_prev  = [eta_hover * ones(8,1); 0; 0];
delta_aero = zeros(3,1);                         % Deflexão de superfícies de controle [rad]


%% CONFIGURAÇÃO DO TEMPO E PRÉ-ALOCAÇÃO DE MEMÓRIA
% Define o vetor de tempo discreto e reserva espaço contíguo na memória (RAM)
% para o histórico de dados (struct 'hist'), evitando lentidão no laço.
t = 0:p.Ts:p.t_sim;
N = length(t);

% --- Estados da Planta e Referências ---
hist.r       = zeros(3, N);     % Posição real [m]
hist.r_bar   = zeros(3, N);     % Posição de referência (Guiamento) [m]
hist.v       = zeros(3, N);     % Velocidade linear [m/s]
hist.q       = zeros(4, N);     % Atitude em quatérnio
hist.q_bar   = zeros(4, N);     % Atitude de referência comandada
hist.w       = zeros(3, N);     % Velocidade angular [rad/s]

% --- Atuadores e Controle ---
hist.varpi   = zeros(p.n_r, N); % Rotação real dos motores [rad/s]
hist.eta     = zeros(p.n_r, N); % Comando enviado aos motores [0 a 1]
hist.f_cmd   = zeros(3, N);     % Força resultante comandada [N]
hist.tau_cmd = zeros(3, N);     % Torque resultante comandado [N.m]
hist.f_star  = zeros(p.n_r, N); % Empuxo ideal exigido por rotor [N]

% --- Diagnóstico Aerodinâmico e Físico ---
hist.alpha   = zeros(1, N);     % Ângulo de ataque [rad]
hist.beta    = zeros(1, N);     % Ângulo de derrapagem [rad]
hist.pdin    = zeros(1, N);     % Pressão dinâmica [Pa]
hist.f_real  = zeros(3, N);     % Força fisicamente realizada [N]
hist.tau_real= zeros(3, N);     % Torque fisicamente realizado [N.m]


%% LAÇO PRINCIPAL DE SIMULAÇÃO (DISCRETE-TIME LOOP)
disp('>> Simulando...');

for k = 1:N
    % 1. Log dos Estados Atuais
    hist.r(:, k)     = mav.r;
    hist.v(:, k)     = mav.v;
    hist.q(:, k)     = mav.q; 
    hist.w(:, k)     = mav.w;
    hist.varpi(:, k) = mav.varpi;
    
    % 2. Atualização de Perturbações e Esforços Externos (Aerodinâmica)
    mav.updateDisturbances(eta_prev, delta_aero);
    
    % 3. Guiamento (Geração de Trajetória)
    ref = guide.getCommand(mav.r);
    hist.r_bar(:, k) = ref.r_bar;
    
    % 4. Controle (Cálculo de Ação e Alocação)
    [eta, log_ctrl] = ctrl.compute(ref, mav);
    
    % 5. Log de Ações do Controlador
    hist.eta(:, k)     = eta; 
    hist.q_bar(:, k)   = D2q(log_ctrl.D_cmd); % Conversão DCM -> Quatérnio
    hist.f_cmd(:, k)   = log_ctrl.f_cmd;
    hist.tau_cmd(:, k) = log_ctrl.tau_cmd;
    hist.f_star(:, k)  = log_ctrl.f_star;
    
    % 6. Planta (Integração RK4 - Propagação da Física)
    mav.integrate(eta, delta_aero);
    eta_prev = eta; % Armazena o comando para causalidade do motor no próximo passo
    
    % --- Extração de Diagnósticos e Pós-Processamento Interno ---
    
    % 7. Diagnóstico Aerodinâmico (Assumindo Vento Inercial = 0)
    D_bg = q2D(mav.q);
    v_body = D_bg * mav.v; 
    V_norm = norm(v_body);
    
    if V_norm > 0.1 % Filtro contra singularidade em baixas velocidades
        hist.alpha(1, k) = atan2(v_body(3), v_body(1));
        hist.beta(1, k)  = asin(max(min(v_body(2)/V_norm, 1), -1));
    else
        hist.alpha(1, k) = 0;
        hist.beta(1, k)  = 0;
    end
    
    rho_ar = 1.225; 
    if isfield(p, 'rho'), rho_ar = p.rho; end 
    hist.pdin(1, k) = 0.5 * rho_ar * V_norm^2;
    
    % 8. Auditoria Física (Esforços de fato realizados pelos atuadores)
    f_rotors_real = p.kf .* (mav.varpi.^2); 
    wrench_real = p.G * f_rotors_real; 
    
    hist.f_real(:, k)   = wrench_real(1:3); 
    hist.tau_real(:, k) = wrench_real(4:6); 
end


%% PÓS-PROCESSAMENTO E EXPORTAÇÃO
% Salva os resultados para análise e aciona o módulo de plotagem
save('Outputs/sim_data.mat', 'hist', 't', 'p');
disp('>> Simulação concluída. Dados salvos em Outputs/sim_data.mat');

% Compila matriz para exportação externa (ex: Python, Excel)
dados_exportacao = [t', hist.r', hist.r_bar', hist.v', hist.q', hist.q_bar', ...
                    hist.w', hist.f_cmd', hist.tau_cmd', hist.f_star', ...
                    hist.varpi', hist.eta'];
                
writematrix(dados_exportacao, 'Outputs/sim_data.csv');

% Geração de Gráficos
plotter.plotAll(t, hist);