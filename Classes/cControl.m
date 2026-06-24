classdef cControl < handle
    properties
        % Físico e Inércia
        m             % massa total da aeronave [kg]
        g             % aceleração da gravidade [m/s^2]
        Jt            % matriz de inércia total (corpo + rotores) [kg.m^2]
        
        % Ganhos de Controle
        K1_pos; K2_pos; % ganhos da malha de posição (Proporcional e Derivativo)
        K1_att; K2_att; % ganhos da malha de atitude (Proporcional e Derivativo)
        
        % Matriz de Alocação
        Theta_pinv    % pseudo-inversa da matriz de alocação de controle (G^+)
        
        f_min; f_max; % limites de empuxo mínimo e máximo exigido por rotor [N]
        kf            % coeficiente de empuxo dos rotores (k_f)
        km            % ganho dos motores (k_m)
        
        % Limites de Segurança do Controlador
        max_tilt_ang  % saturação de Pitch/Roll em queda livre [rad]
        Tz_min        % limite para denominador em atan2 [N]
        trace_min     % limite inferior de traço para vetor Gibbs
    end
    
    methods   
        %% CONSTRUTORA
        % Inicializa a controladora carregando os parâmetros físicos, ganhos 
        % de malha e calculando a matriz pseudo-inversa de alocação de controle.
        function obj = cControl(sControl)
            % Parâmetros Físicos e Inerciais
            obj.m  = sControl.m; 
            obj.g  = sControl.g; 
            obj.Jt = sControl.Jt;
            
            % Ganhos das Malhas de Controle
            obj.K1_pos = sControl.K1_pos; 
            obj.K2_pos = sControl.K2_pos;
            obj.K1_att = sControl.K1_att; 
            obj.K2_att = sControl.K2_att;
            
            % Parâmetros e Limites dos Atuadores
            obj.f_min = sControl.f_min; 
            obj.f_max = sControl.f_max; 
            obj.kf    = sControl.kf; 
            obj.km    = sControl.km;
            
            % Limites de Segurança do Controlador
            obj.max_tilt_ang = sControl.max_tilt_ang;
            obj.Tz_min       = sControl.Tz_min;
            obj.trace_min    = sControl.trace_min;
            
            % Pseudo-inversa calculada dinamicamente com base na matriz G configurada
            obj.Theta_pinv = pinv(sControl.G); 
        end
        
        
        %% COMPUTAÇÃO DA LEI DE CONTROLE
        % Executa as malhas de controle em cascata (Posição -> Atitude),
        % resolve a alocação híbrida e converte os esforços em comandos 
        % de atuação (eta) para os motores. O método é desacoplado da planta.
        function [eta, log_ctrl] = compute(obj, ref, r, v, q, w, f_aero, tau_aero)
            % --- Desarmamento e Corte de Motores no Pouso ---
            if ref.flight_mode == 3
                eta = zeros(10, 1);
                log_ctrl.D_cmd   = eye(3);
                log_ctrl.f_cmd   = zeros(3,1);
                log_ctrl.tau_cmd = zeros(3,1);
                log_ctrl.f_star  = zeros(10,1);
                return;
            end
            
            % Extração da Matriz de Atitude Atual (quaternião com escalar no final)
            D_bg = q2D(q);
            
            % --- Malha de Controle de Posição ---
            
            % 1. Aceleração de Comando e Força Inercial Resultante
            e_r = r - ref.r_bar;
            e_v = v - ref.v_bar;
            
            a_com = ref.a_bar - obj.K1_pos * e_r - obj.K2_pos * e_v;
            f_I_cmd = obj.m * (a_com + [0; 0; obj.g]) - D_bg' * f_aero;
            
            % --- Mapeamento Híbrido (CCA - Veículo Completamente Atuado) ---
            
            % 2. Decomposição de Esforços no Referencial de Guinada (Psi)
            % Correção: Como rotz() é uma matriz passiva, ela já mapeia do Inercial para o Corpo. Não usar transposta!
            R_z = rotz(ref.alpha_bar(3)); 
            f_psi = R_z * f_I_cmd;  
            
            Tx_req = f_psi(1); 
            Ty_req = f_psi(2);
            Tz_req = f_psi(3); 
            
            % 3. Lógica Híbrida e Desacoplamento Geométrico Exato
            max_ang = obj.max_tilt_ang;
            if ref.flight_mode == 1 || Tx_req < 0
                % FREAR ou MULTICÓPTERO: O drone usa rotores verticais para tudo (Pitch UP + Roll)
                Tx_alloc = 0;
                Tz_alloc = norm([Tx_req, Ty_req, Tz_req]); % Empuxo vertical total real
                
                theta_des = max(min(atan2(Tx_req, max(Tz_alloc, obj.Tz_min)), max_ang), -max_ang);
                phi_des   = max(min(atan2(-Ty_req, max(Tz_req, obj.Tz_min)), max_ang), -max_ang);
            else
                % ACELERAR (Híbrido): Rotores frontais lidam com Tx. Verticais lidam com Ty e Tz.
                Tx_alloc = Tx_req;
                Tz_alloc = norm([0, Ty_req, Tz_req]); % Empuxo vertical total real
                
                theta_des = 0;
                phi_des   = max(min(atan2(-Ty_req, max(Tz_req, obj.Tz_min)), max_ang), -max_ang);
            end
            
            % 4. Matriz de Comando Passiva Desejada (Solo -> Corpo)
            D_cmd = rotx(phi_des) * roty(theta_des) * rotz(ref.alpha_bar(3));
            
            % --- Malha de Controle de Atitude ---
            
            % 5. Erro de Atitude (Mapeado no Vetor de Gibbs para SO(3))
            D_tilde = D_bg * D_cmd';
            den = max(1 + trace(D_tilde), obj.trace_min);
            e_alpha = (1 / den) * [D_tilde(2,3) - D_tilde(3,2);
                                   D_tilde(3,1) - D_tilde(1,3);
                                   D_tilde(1,2) - D_tilde(2,1)];
            
            % 6. Controle PD Geométrico e Cancelamento Giroscópico de Rotação
            alpha_com = - obj.K1_att * e_alpha - obj.K2_att * (w - ref.w_bar);
            tau_b_c = obj.Jt * alpha_com + skew(w) * (obj.Jt * w) - tau_aero;
            
            % --- Alocação de Controle e Atuadores ---
            
            % 7. Aplicação da Pseudo-Inversa e Saturação Física dos Rotores
            u = [Tx_alloc; 0; Tz_alloc; tau_b_c]; 
            f_irrestrito = obj.Theta_pinv * u;
            
            f_star = min(max(f_irrestrito, obj.f_min), obj.f_max);
            eta = sqrt(f_star ./ (obj.kf .* (obj.km.^2)));
            
            % --- Encapsulamento ---
            
            % 8. Pacote de Logs de Controle (Para Análise e Plotagem)
            log_ctrl.D_cmd   = D_cmd;                   % Matriz de atitude comandada
            log_ctrl.f_cmd   = [Tx_alloc; 0; Tz_alloc]; % Força comandada no corpo (3x1) [N]
            log_ctrl.tau_cmd = tau_b_c;                 % Torque comandado (3x1) [N.m]
            log_ctrl.f_star  = f_star;                  % Empuxo individual exigido após alocação [N]
        end
        
    end
end