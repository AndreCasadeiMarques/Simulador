%% =====================================================================
% PARÂMETROS DE SIMULAÇÃO - eVTOL GD-350
% Descrição: Define todos os parâmetros físicos, aerodinâmicos, 
% geométricos e de controle organizados em structs isoladas por domínio.
% =====================================================================

%% 1. STRUCTS AUXILIARES E PARÂMETROS GERAIS
m_val       = 390.0;                      % massa total da aeronave [kg]
g_val       = 9.80665;                    % aceleração da gravidade [m/s^2]
n_r_val     = 10;                         % número total de rotores
Ts_val      = 0.001;                      % passo de integração [s]
t_sim_val   = 190.0;                       % tempo total de simulação [s]

% Matriz de inércia do corpo (sem os rotores) [kg.m^2]
Jb = diag([258.9, 331.8, 575.4]); 
Jb(1,3) = -14.0; 
Jb(3,1) = -14.0;

% Parâmetros dos Motores e Rotores
Jr_val = 0.12 * ones(n_r_val, 1);
mum_val = 0.25 * ones(n_r_val, 1);
km_val = 1000 * ones(n_r_val, 1);
kf_val = (1115 / 1000^2) * ones(n_r_val, 1);
k_val = 0.0683 * ones(n_r_val, 1);

w_min_val = [0 * ones(8, 1); 0; 0];
w_max_val = 1000 * ones(n_r_val, 1);
f_min_val = [0 * ones(8, 1); 0; 0];
f_max_val = 1115 * ones(n_r_val, 1);

% Geometria dos Rotores (Posição em relação ao CG)
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

sigma_val = [1; 1; 1; 1; -1; -1; -1; -1; 1; -1];
gamma_v = [-3, -3, -3, -3, 3, 3, 3, 3] * (pi/180);

% Pré-computação das matrizes geométricas dos rotores
G_val = zeros(6, n_r_val); 
D_rb_val = zeros(3, 3, n_r_val); 
Js_val = zeros(3, 3, n_r_val);

for i = 1:8
    x_i  = ell_b(1, i); 
    y_i  = ell_b(2, i); 
    d_xy = sqrt(x_i^2 + y_i^2); 
    th   = gamma_v(i);
    
    D_rib_T = [ -(x_i/d_xy)*cos(th),  (y_i/d_xy), -(x_i/d_xy)*sin(th);
                -(y_i/d_xy)*cos(th), -(x_i/d_xy), -(y_i/d_xy)*sin(th);
                -sin(th),             0,           cos(th) ];
                
    D_rb_val(:,:,i) = D_rib_T;
    Js_val(:,:,i)   = diag([0.05, 0.05, Jr_val(i)]); 
    
    gamma_f = D_rib_T * [0; 0; 1];
    G_val(:, i) = [gamma_f; cross(ell_b(:,i), gamma_f) - k_val(i)*sigma_val(i)*gamma_f];
end

for i = 9:10
    D_rib_T = [ 0, 0, 1;
                0, 1, 0;
               -1, 0, 0];
              
    D_rb_val(:,:,i) = D_rib_T;
    Js_val(:,:,i)   = diag([0.05, 0.05, Jr_val(i)]);
    
    gamma_f = [1; 0; 0]; 
    G_val(:, i) = [gamma_f; cross(ell_b(:,i), gamma_f) - k_val(i)*sigma_val(i)*gamma_f];
end

% Composição do Tensor de Inércia Total (Corpo + Rotores)
Jb_rot = zeros(3,3);
for i = 1:n_r_val
    Jb_rot = Jb_rot + D_rb_val(:,:,i) * Js_val(:,:,i) * D_rb_val(:,:,i)'; 
end
Jt_val = Jb + Jb_rot; 
Jt_inv_val = inv(Jt_val);


%% 2. CONFIGURAÇÃO DA STRUCT S_MAV (PLANTA)
sMav = struct();
sMav.m      = m_val;
sMav.g      = g_val;
sMav.Jt     = Jt_val;
sMav.Jt_inv = Jt_inv_val;
sMav.Ts     = Ts_val;

sMav.mum    = mum_val;
sMav.km     = km_val;
sMav.kf     = kf_val;
sMav.w_min  = w_min_val;
sMav.w_max  = w_max_val;
sMav.n_r    = n_r_val;
sMav.sigma  = sigma_val;
sMav.D_rb   = D_rb_val;
sMav.Js     = Js_val;
sMav.G      = G_val;

% Parâmetros Aerodinâmica
sMav.rho = 1.225; % densidade do ar [kg/m^3]
sMav.Aa  = 2.0;   % área de referência [m^2]
sMav.c   = 0.5;   % corda média aerodinâmica [m]

% Limites e Parâmetros de Diagnóstico/Física
sMav.v_aero_min = 2.0; % velocidade mínima para cálculo aerodinâmico [m/s]
sMav.v_diag_min = 0.1; % velocidade mínima para diagnóstico de ângulos [m/s]
sMav.alpha_max  = 20 * (pi/180); % limite de saturação numérico de ataque [rad]
sMav.beta_max   = 10 * (pi/180); % limite de saturação numérico de derrapagem [rad]
sMav.ground_friction = 0.9;      % fator de decaimento de velocidade XY no solo
sMav.gamma_v    = 3 * (pi/180);  % cant angle dos rotores verticais [rad]

sMav.CD0 =  0.0312;   sMav.CDa =  0.0;      sMav.CDq = -0.5926;   sMav.CDde =  0.0084;
sMav.CYb = -0.4727;   sMav.CYp =  0.0958;   sMav.CYr =  0.1665;   sMav.CYda =  0.0;      sMav.CYdr =  0.0034;
sMav.CL0 =  0.0;      sMav.CLa =  5.8392;   sMav.CLq = 10.2236;   sMav.CLde =  0.0084;
sMav.Clb = -0.0312;   sMav.Clp = -0.5926;   sMav.Clr =  0.2390;   sMav.Clda =  0.0045;   sMav.Cldr = -0.0;
sMav.Cm0 =  0.0;      sMav.Cma = -1.7199;   sMav.Cmq = -21.9187;  sMav.Cmde = -0.0309;
sMav.Cnb =  0.0726;   sMav.Cnp = -0.0810;   sMav.Cnr = -0.0732;   sMav.Cnda =  0.0001;   sMav.Cndr = -0.0010;


%% 3. CONFIGURAÇÃO DA STRUCT S_CONTROL (CONTROLADOR)
sControl = struct();
sControl.m      = m_val;
sControl.g      = g_val;
sControl.Jt     = Jt_val;
sControl.f_min  = f_min_val;
sControl.f_max  = f_max_val;
sControl.kf     = kf_val;
sControl.km     = km_val;
sControl.G      = G_val;

% Parâmetros e Limites de Segurança do Controlador
sControl.max_tilt_ang = 30 * (pi/180); % saturação de Pitch/Roll em queda livre [rad]
sControl.Tz_min       = 0.1;           % limite para denominador em atan2 [N]
sControl.trace_min    = 1e-4;          % limite inferior de traço para vetor Gibbs
% Ganhos de Controle e Sintonia (wn e zeta)
wn_pos = 0.2;
wn_att = 1.0;
zeta   = 1.0;

sControl.K1_pos = diag([wn_pos^2, wn_pos^2, wn_pos^2]);
sControl.K2_pos = diag([2*zeta*wn_pos, 2*zeta*wn_pos, 2*zeta*wn_pos]);
sControl.K1_att = diag([wn_att^2, wn_att^2, wn_att^2]);
sControl.K2_att = diag([2*zeta*wn_att, 2*zeta*wn_att, 2*zeta*wn_att]);


%% 4. CONFIGURAÇÃO DA STRUCT S_GUIDANCE (GUIAMENTO)
sGuidance = struct();
sGuidance.Ts      = Ts_val;
sGuidance.mode    = 'multicoptero'; % 'armado', 'multicoptero', 'transicao'
sGuidance.W_r     = [0,0,0; 0,0,100; 100,0,100; 100,100,100; 0,100,100; 0,0,100; 0,0,0]';
sGuidance.W_alpha = zeros(3, size(sGuidance.W_r, 2));
sGuidance.v_avg   = 5.0; % velocidade média do Minimum Jerk de decolagem [m/s]
sGuidance.v_avg_landing = 5.0; % velocidade média do Minimum Jerk de pouso [m/s]
sGuidance.R_acc   = 5.0; % raio de wayset de cruzeiro [m]
sGuidance.R_acc_landing = 0.8; % raio de wayset apertado para o pré-pouso [m]
sGuidance.v_max   = 5.0; % velocidade máxima de cruzeiro [m/s]
sGuidance.a_max   = 0.5; % aceleração máxima [m/s^2] (Reduzido para o drone real acompanhar)
sGuidance.wn_ref  = 0.2; % frequência do campo vetorial (Reduzido para parear com wn_pos)

% Limites de Segurança e Tolerância do Guiamento
sGuidance.v_stop_tol   = 0.5; % velocidade para considerar parada total [m/s]
sGuidance.yaw_rate_max = 30 * (pi/180); % saturação de velocidade de guinada [rad/s]
sGuidance.tj_min       = 0.1; % tempo mínimo de trajetória Minimum Jerk [s]


%% 5. SALVAMENTO E EXPORTAÇÃO
sSim.t_sim = t_sim_val;
sSim.Ts = Ts_val;
sSim.n_r = n_r_val;
save('parameters.mat', 'sMav', 'sControl', 'sGuidance', 'sSim');