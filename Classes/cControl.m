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
        
        % Limites e Parâmetros dos Atuadores
        f_min; f_max; % limites de empuxo mínimo e máximo exigido por rotor [N]
        kf            % coeficiente de empuxo dos rotores (k_f)
        km            % ganho dos motores (k_m)
    end
    
    methods   
        %% CONSTRUTORA
        % Inicializa a controladora carregando os parâmetros físicos, ganhos 
        % de malha e calculando a matriz pseudo-inversa de alocação de controle.
        function obj = cControl(p)
            % Parâmetros Físicos e Inerciais
            obj.m  = p.m; 
            obj.g  = p.g; 
            obj.Jt = p.Jt;
            
            % Ganhos das Malhas de Controle
            obj.K1_pos = p.K1_pos; 
            obj.K2_pos = p.K2_pos;
            obj.K1_att = p.K1_att; 
            obj.K2_att = p.K2_att;
            
            % Parâmetros e Limites dos Atuadores
            obj.f_min = p.f_min; 
            obj.f_max = p.f_max; 
            obj.kf    = p.kf; 
            obj.km    = p.km;
            
            % Alocador R^5 (Ignora a força em Y, que é resolvida por Rolagem)
            % Nota: A pseudoinversa distribui os esforços perfeitamente na planta sobreatuada
            obj.Theta_pinv = pinv(p.G); 
        end
        
        
        %% COMPUTAÇÃO DA LEI DE CONTROLE
        % Executa as malhas de controle em cascata (Posição -> Atitude), 
        % resolve a alocação híbrida (transição avião/helicóptero) e 
        % converte os esforços em comandos de atuação (eta) para os motores.
        function [eta, log_ctrl] = compute(obj, ref, mav)
            % Extração da Matriz de Atitude Atual
            D_bg = q2D(mav.q);
            
            % --- Malha de Controle de Posição ---
            
            % 1. Aceleração de Comando e Força Inercial Resultante
            e_r = mav.r - ref.r_bar;
            e_v = mav.v - ref.v_bar;
            
            a_com = ref.a_bar - obj.K1_pos * e_r - obj.K2_pos * e_v;
            f_I_cmd = obj.m * (a_com + [0; 0; obj.g]) - D_bg' * mav.f_aero;
            
            % --- Mapeamento Híbrido (CCA - Veículo Completamente Atuado) ---
            
            % 2. Decomposição de Esforços no Referencial de Guinada (Psi)
            R_z = rotz(ref.alpha_bar(3)); 
            f_psi = R_z * f_I_cmd;  
            
            Tx_req = f_psi(1); 
            Ty_req = f_psi(2);
            Tz_req = f_psi(3); 
            
            % 3. Lógica Híbrida Longitudinal (Pitch e Thrust Frontal)
            if Tx_req < 0
                % FREAR: O drone empina (Nose UP) e corta os rotores horizontais
                theta_des = asin(max(min(Tx_req / max(Tz_req, 0.1), 1), -1));
                Tx_alloc  = 0; 
            else
                % ACELERAR: O drone voa nivelado e usa os rotores horizontais
                theta_des = 0; 
                Tx_alloc  = Tx_req; 
            end
            
            % 4. Esforço Lateral (Sempre resolvido mecanicamente por Roll LEFT/RIGHT)
            phi_des = asin(max(min(-Ty_req / max(Tz_req, 0.1), 1), -1));
            
            % Matriz de Comando Passiva Desejada (Solo -> Corpo)
            D_cmd = rotx(phi_des) * roty(theta_des) * rotz(ref.alpha_bar(3));
            
            % --- Malha de Controle de Atitude ---
            
            % 5. Erro de Atitude (Mapeado no Vetor de Gibbs para SO(3))
            D_tilde = D_bg * D_cmd';
            den = max(1 + trace(D_tilde), 1e-4);
            e_alpha = (1 / den) * [D_tilde(2,3) - D_tilde(3,2);
                                   D_tilde(3,1) - D_tilde(1,3);
                                   D_tilde(1,2) - D_tilde(2,1)];
            
            % 6. Controle PD Geométrico e Cancelamento Giroscópico de Rotação
            alpha_com = - obj.K1_att * e_alpha - obj.K2_att * (mav.w - ref.w_bar);
            tau_b_c = obj.Jt * alpha_com + skew(mav.w) * (obj.Jt * mav.w) - mav.tau_aero;
            
            % --- Alocação de Controle e Atuadores ---
            
            % 7. Aplicação da Pseudo-Inversa e Saturação Física dos Rotores
            u = [Tx_alloc; 0; Tz_req; tau_b_c]; 
            f_irrestrito = obj.Theta_pinv * u;
            
            f_star = min(max(f_irrestrito, obj.f_min), obj.f_max);
            eta = sqrt(f_star ./ (obj.kf .* (obj.km.^2)));
            
            % --- Encapsulamento ---
            
            % 8. Pacote de Logs de Controle (Para Análise e Plotagem)
            log_ctrl.D_cmd   = D_cmd;                   % Matriz de atitude comandada
            log_ctrl.f_cmd   = [Tx_alloc; 0; Tz_req];   % Força comandada no corpo (3x1) [N]
            log_ctrl.tau_cmd = tau_b_c;                 % Torque comandado (3x1) [N.m]
            log_ctrl.f_star  = f_star;                  % Empuxo individual exigido após alocação (10x1) [N]
        end
        
    end
end