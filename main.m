% main.m
% Descrição: Laço principal de simulação GNC usando RK4.

clear; clc; close all;

% Adiciona todas as subpastas ao Path do MATLAB para enxergar as Classes e Funções
addpath(genpath(pwd));

% Roda a carrega o script de parâmetros
run('parameters.m');
load('parameters.mat', 'p');

% Instanciação de Objetos (Injeção de Dependência)
mav   = cMav(p);
ctrl  = cControl(p);
guide = cGuidance(p);
plotter = cPlotter();

% 2. Configuração do Laço Temporal
t = 0:p.Ts:p.t_sim;
N = length(t);

% Estrutura de Log
hist.r = zeros(3, N); hist.r_bar = zeros(3, N);
hist.v = zeros(3, N); 
hist.q = zeros(4, N); hist.w = zeros(3, N); 
hist.varpi = zeros(p.n_r, N); hist.eta = zeros(p.n_r, N);

delta_aero = zeros(3,1); 
eta_prev   = zeros(p.n_r, 1);

disp('>> Simulando...');

for k = 1:N
    % Log de Estados Atuais
    hist.r(:, k) = mav.r;
    hist.v(:, k) = mav.v;
    hist.q(:, k) = mav.q; 
    hist.w(:, k) = mav.w;
    hist.varpi(:, k) = mav.varpi;

    % Perturbações
    mav.updateDisturbances(eta_prev, delta_aero);

    % Guiamento
    ref = guide.getCommand(mav.r);
    hist.r_bar(:, k) = ref.r_bar;

    % Controle
    [eta, log_ctrl] = ctrl.compute(ref, mav);
    hist.eta(:, k) = eta; 

    hist.q_bar(:, k)   = D2q(log_ctrl.D_cmd); % Converte DCM para quatérnio
    hist.f_cmd(:, k)   = log_ctrl.f_cmd;
    hist.tau_cmd(:, k) = log_ctrl.tau_cmd;
    hist.f_star(:, k)  = log_ctrl.f_star;

    % Planta (Integração RK4)
    mav.integrate(eta, delta_aero);

    eta_prev = eta; % Causalidade do motor
end

% 4. Pós-Processamento
save('Outputs/sim_data.mat', 'hist', 't', 'p');
disp('>> Simulação concluída. Dados salvos em Outputs/sim_data.mat');
dados_exportacao = [t', hist.r', hist.r_bar', hist.v', hist.q', hist.q_bar', hist.w', hist.f_cmd', hist.tau_cmd', hist.f_star', hist.varpi', hist.eta'];
writematrix(dados_exportacao, 'Outputs/sim_data.csv');
plotter.plotAll(t, hist);