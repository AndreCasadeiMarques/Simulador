%% =====================================================================
% PARÂMETROS DE SIMULAÇÃO - eVTOL GD-350
% Descrição: Define todos os parâmetros físicos, aerodinâmicos, 
% geométricos e de controle para a simulação do veículo.
% =====================================================================
p = struct();

%% SIMULAÇÃO E TEMPO
p.t_sim = 40.0;                         % tempo total de simulação [s]
p.Ts    = 0.002;                        % passo de integração / tempo de amostragem [s]


%% PARÂMETROS FÍSICOS E INERCIAIS
p.m       = 390.0;                      % massa total da aeronave [kg]
p.g       = 9.80665;                    % aceleração da gravidade [m/s^2]

% Matriz de inércia do corpo (sem os rotores) [kg.m^2]
p.Jb      = diag([258.9, 331.8, 575.4]); 
p.Jb(1,3) = -14.0; 
p.Jb(3,1) = -14.0;

% Matriz de inércia base/referência (mantida por legado/verificação) [kg.m^2]
p.Jbatata = [247.8,     0,     0;
                 0, 303.8,     0;
                 0,     0, 540.8];


%% SISTEMA DE PROPULSÃO E ATUADORES
p.n_r   = 10;                                   % número total de rotores
p.Jr    = 0.12 * ones(p.n_r, 1);                % momento de inércia individual do rotor [kg.m^2]
p.mum   = 0.25 * ones(p.n_r, 1);                % constante de tempo mecânica dos motores [s]
p.km    = 1000 * ones(p.n_r, 1);                % ganho estático dos motores (comando para RPM)

% Coeficientes de Esforço
p.kf    = (1114.46 / 1000^2) * ones(p.n_r, 1);  % coeficiente de empuxo dos rotores [N/(rad/s)^2]
p.k     = 0.02 * ones(p.n_r, 1);                % razão torque/empuxo (k_tau/k_f) [m]

% Limites de Saturação Física
p.w_min = [100 * ones(8, 1); 0; 0];             % velocidade angular mínima [rad/s]
p.w_max = 1000 * ones(p.n_r, 1);                % velocidade angular máxima [rad/s]
p.f_min = [10 * ones(8, 1); 0; 0];              % empuxo mínimo exigido por rotor [N]
p.f_max = 1114.46 * ones(p.n_r, 1);             % empuxo máximo gerado por rotor [N]


%% GEOMETRIA DOS ROTORES
% Posição vetorial [x; y; z] de cada rotor em relação ao CG da aeronave [m]
ell_b = [ 1.531, -1.305, 0.506;  
          1.531,  1.305, 0.506; 
         -1.531,  1.305, 0.506; 
         -1.531, -1.305, 0.506;
          1.531,  1.305, 0.244;  
          1.531, -1.305, 0.244; 
         -1.531, -1.305, 0.244; 
         -1.531,  1.305, 0.244;
         -1.556,  0.000, 0.366;  
          1.556,  0.000, 0.366]';

% Sentido de rotação (+1 anti-horário, -1 horário)
p.sigma = [1; 1; 1; 1; -1; -1; -1; -1; 1; -1];

% Ângulo de diedro (inclinação lateral) dos rotores de voo vertical (1 a 8) [rad]
gamma_v = [-3, -3, -3, -3, 3, 3, 3, 3] * (pi/180);


%% PRÉ-COMPUTAÇÃO MATRICIAL (ALOCAÇÃO E INÉRCIA TOTAL)
% Inicialização das matrizes tridimensionais e do alocador
p.G    = zeros(6, p.n_r); 
p.D_rb = zeros(3, 3, p.n_r); 
p.Js   = zeros(3, 3, p.n_r);

% 1. Configuração dos Rotores de Sustentação Vertical (1 a 8) com Diedro
for i = 1:8
    x_i  = ell_b(1, i); 
    y_i  = ell_b(2, i); 
    d_xy = sqrt(x_i^2 + y_i^2); 
    th   = gamma_v(i);
    
    % Matriz analítica exata de atitude local do rotor (Transposta D^{r_i/b}^T)
    D_rib_T = [ -(x_i/d_xy)*cos(th),  (y_i/d_xy), -(x_i/d_xy)*sin(th);
                -(y_i/d_xy)*cos(th), -(x_i/d_xy), -(y_i/d_xy)*sin(th);
                -sin(th),             0,           cos(th) ];
                
    p.D_rb(:,:,i) = D_rib_T;
    p.Js(:,:,i)   = diag([0.05, 0.05, p.Jr(i)]); 
    
    % Direção do vetor de empuxo (3ª coluna exata da matriz acima)
    gamma_f = D_rib_T * [0; 0; 1];
    
    % Montagem do Alocador G (com correção do torque reverso de guinada)
    p.G(:, i) = [gamma_f; cross(ell_b(:,i), gamma_f) - p.k(i)*p.sigma(i)*gamma_f];
end

% 2. Configuração dos Rotores Horizontais (9 e 10) para Empuxo Frontal
for i = 9:10
    D_rib_T = [ 0, 0, 1;
                0, 1, 0;
               -1, 0, 0];
              
    p.D_rb(:,:,i) = D_rib_T;
    p.Js(:,:,i)   = diag([0.05, 0.05, p.Jr(i)]);
    
    % Motores horizontais apontam estritamente para o eixo X
    gamma_f = [1; 0; 0]; 
    
    p.G(:, i) = [gamma_f; cross(ell_b(:,i), gamma_f) - p.k(i)*p.sigma(i)*gamma_f];
end

% 3. Composição do Tensor de Inércia Total (Corpo + Rotores)
Jb_rot = zeros(3,3);
for i = 1:p.n_r
    Jb_rot = Jb_rot + p.D_rb(:,:,i) * p.Js(:,:,i) * p.D_rb(:,:,i)'; 
end

p.Jt     = p.Jb + Jb_rot; 
p.Jt_inv = inv(p.Jt);


%% PARÂMETROS AERODINÂMICOS
p.rho = 1.225;          % densidade do ar [kg/m^3]
p.Aa  = 7.34;           % área de asa/referência aerodinâmica [m^2]
p.c   = 0.7917;         % corda de referência [m]

% Coeficientes de Arrasto (Eixo X)
p.CD0 =  0.0312;   p.CDa =  0.0;      p.CDq = -0.5926;   p.CDde =  0.0084;

% Coeficientes de Força Lateral (Eixo Y)
p.CYb = -0.4727;   p.CYp =  0.0958;   p.CYr =  0.1665;   p.CYda =  0.0;      p.CYdr =  0.0034;

% Coeficientes de Sustentação (Eixo Z)
p.CL0 =  0.0;      p.CLa =  5.8392;   p.CLq = 10.2236;   p.CLde =  0.0084;

% Coeficientes de Momento de Rolagem (Roll)
p.Clb = -0.0312;   p.Clp = -0.5926;   p.Clr =  0.2390;   p.Clda =  0.0045;   p.Cldr = -0.0;

% Coeficientes de Momento de Arfagem (Pitch)
p.Cm0 =  0.0;      p.Cma = -1.7199;   p.Cmq = -21.9187;  p.Cmde = -0.0309;

% Coeficientes de Momento de Guinada (Yaw)
p.Cnb =  0.0726;   p.Cnp = -0.0810;   p.Cnr = -0.0732;   p.Cnda =  0.0001;   p.Cndr = -0.0010;


%% GANHOS DE CONTROLE E SINTONIA
wn_pos = 0.1;           % frequência natural da malha de posição [rad/s]
wn_att = 0.5;           % frequência natural da malha de atitude [rad/s]
zeta   = 1.0;           % fator de amortecimento (crítico) para as malhas

p.K1_pos = diag([wn_pos^2, wn_pos^2, wn_pos^2]);
p.K2_pos = diag([2*zeta*wn_pos, 2*zeta*wn_pos, 2*zeta*wn_pos]);

p.K1_att = diag([wn_att^2, wn_att^2, wn_att^2]);
p.K2_att = diag([2*zeta*wn_att, 2*zeta*wn_att, 2*zeta*wn_att]);


%% GUIAMENTO (DRONE VIRTUAL) E MISSÃO
% Matrizes de Waypoints da Missão
p.W_r     = [0,0,0; 0,0,10; 100,0,20]'; % waypoints de posição 3xN [m]
p.W_alpha = zeros(3, 3);                % waypoints de atitude Euler 3xN [rad]

% Parâmetros do Gerador de Trajetória Fantasma
p.R_acc    = 0.5;       % raio esférico de aceitação do waypoint [m]
p.v_max    = 10.0;      % saturação de velocidade do drone fantasma [m/s]
p.a_max    = 2.0;       % saturação de aceleração do drone fantasma [m/s^2]
p.wn_ref   = 0.8;       % frequência natural de rastreio da malha virtual [rad/s]
p.zeta_ref = 1.0;       % amortecimento da malha virtual


%% SALVAMENTO
save('parameters.mat', 'p');