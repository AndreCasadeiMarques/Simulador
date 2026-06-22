% parameters.m
p = struct();

% Físicos
p.m   = 390.0;          
p.g   = 9.80665;        
p.Ts  = 0.002;          

p.Jbatata = [247.8      0      0;
            0  303.8      0;
            0      0  540.8];
p.Jb  = diag([258.9, 331.8, 575.4]); 
p.Jb(1,3) = -14.0; p.Jb(3,1) = -14.0;

% Atuadores (Limites movidos para cá - Fim do Hardcode)
p.n_r   = 10;                     
p.mum   = 0.25 * ones(p.n_r, 1);  
p.km    = 1000 * ones(p.n_r, 1);  
p.w_min = [100 * ones(8, 1); 0; 0];                    
p.w_max = 1000 * ones(p.n_r, 1);                   
p.f_min = [10 * ones(8, 1); 0; 0];      % Empuxo mínimo [N]
p.f_max = 1114.46 * ones(p.n_r, 1);     % Empuxo máximo [N]
p.kf    = (1114.46 / 1000^2) * ones(p.n_r, 1); 
p.k     = 0.02 * ones(p.n_r, 1);  
p.Jr    = 0.12 * ones(p.n_r, 1); 

% Geometria
ell_b = [1.531, -1.305, 0.506; 1.531, 1.305, 0.506; -1.531, 1.305, 0.506; -1.531, -1.305, 0.506;
    1.531, 1.305, 0.244; 1.531, -1.305, 0.244; -1.531, -1.305, 0.244; -1.531, 1.305, 0.244;
    -1.556, 0.000, 0.366; 1.556, 0.000, 0.366]';
p.sigma = [1; 1; 1; 1; -1; -1; -1; -1; 1; -1];

% Matrizes de Alocação e Tensores
p.G = zeros(6, p.n_r); p.D_rb = zeros(3, 3, p.n_r); p.Js = zeros(3, 3, p.n_r);
gamma_v = [-3, -3, -3, -3, 3, 3, 3, 3] * (pi/180);

for i = 1:8
    x_i = ell_b(1, i); y_i = ell_b(2, i); d_xy = sqrt(x_i^2 + y_i^2); th = gamma_v(i);

    % Matriz Analítica EXATA do seu relatório (Transposta D^{r_i/b}^T)
    D_rib_T = [ -(x_i/d_xy)*cos(th),  (y_i/d_xy), -(x_i/d_xy)*sin(th);
                -(y_i/d_xy)*cos(th), -(x_i/d_xy), -(y_i/d_xy)*sin(th);
                -sin(th),             0,           cos(th) ];

    p.D_rb(:,:,i) = D_rib_T;
    p.Js(:,:,i)   = diag([0.05, 0.05, p.Jr(i)]); 

    % Gamma_f é a 3ª coluna exata da matriz acima (conforme o relatório)
    gamma_f = D_rib_T * [0; 0; 1];

    % Montagem do Alocador G (com a correção do torque reverso de guinada)
    p.G(:, i) = [gamma_f; cross(ell_b(:,i), gamma_f) - p.k(i)*p.sigma(i)*gamma_f];
end

for i = 9:10
    D_rib_T = [0, 0, 1;
               0, 1, 0;
              -1, 0, 0];
              
    p.D_rb(:,:,i) = D_rib_T;
    p.Js(:,:,i)   = diag([0.05, 0.05, p.Jr(i)]);
    
    gamma_f = [1; 0; 0]; % Motores horizontais apontam para o eixo X
    
    p.G(:, i) = [gamma_f; cross(ell_b(:,i), gamma_f) - p.k(i)*p.sigma(i)*gamma_f];
end

Jb_rot = zeros(3,3);
for i = 1:p.n_r
    Jb_rot = Jb_rot + p.D_rb(:,:,i) * p.Js(:,:,i) * p.D_rb(:,:,i)'; 
end
p.Jt = p.Jb + Jb_rot; 
p.Jt_inv = inv(p.Jt);

% Aerodinâmica
p.rho = 1.225; 
p.Aa = 7.34; 
p.c = 0.7917;         
p.CD0 = 0.0312; p.CDa = 0.0; p.CDq = -0.5926; p.CDde = 0.0084;
p.CYb = -0.4727; p.CYp = 0.0958; p.CYr = 0.1665; p.CYda = 0.0; p.CYdr = 0.0034;
p.CL0 = 0.0; p.CLa = 5.8392; p.CLq = 10.2236; p.CLde = 0.0084;
p.Clb = -0.0312; p.Clp = -0.5926; p.Clr = 0.2390; p.Clda = 0.0045; p.Cldr = -0.0;
p.Cm0 = 0.0; p.Cma = -1.7199; p.Cmq = -21.9187; p.Cmde = -0.0309;
p.Cnb = 0.0726; p.Cnp = -0.0810; p.Cnr = -0.0732; p.Cnda = 0.0001; p.Cndr = -0.0010;

% Ganhos
wn_pos = 0.1; wn_att = 0.5; zeta = 1.0; 
p.K1_pos = diag([wn_pos^2, wn_pos^2, wn_pos^2]);
p.K2_pos = diag([2*zeta*wn_pos, 2*zeta*wn_pos, 2*zeta*wn_pos]);
p.K1_att = diag([wn_att^2, wn_att^2, wn_att^2]);
p.K2_att = diag([2*zeta*wn_att, 2*zeta*wn_att, 2*zeta*wn_att]);

% Simulação e Guiamento
p.t_sim = 40.0;
p.W_r = [0,0,0; 0,0,10; 100,0,20]';
p.W_alpha = zeros(3, 3);
p.R_acc = 0.5; p.v_max = 10.0; p.a_max = 2.0; 
p.wn_ref = 0.8; p.zeta_ref = 1.0;

save('parameters.mat', 'p');