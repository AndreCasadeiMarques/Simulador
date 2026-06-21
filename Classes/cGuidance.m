classdef cGuidance < handle
    properties
        W_r, W_alpha, R_acc, v_max, a_max, wn_ref, Ts
        idx, n_w, r_ref, v_ref, a_ref
        alpha_ref, w_ref
    end
    
    methods
        function obj = cGuidance(p)
            obj.W_r = p.W_r; 
            obj.W_alpha = p.W_alpha;
            obj.R_acc = p.R_acc; 
            obj.v_max = p.v_max; 
            obj.a_max = p.a_max;
            obj.wn_ref = p.wn_ref; 
            obj.Ts = p.Ts;
            
            obj.idx = 1; 
            obj.n_w = size(obj.W_r, 2);
            
            % O drone virtual nasce exatamente no primeiro waypoint
            obj.r_ref = obj.W_r(:, 1); 
            obj.v_ref = zeros(3,1); 
            obj.a_ref = zeros(3,1);
            
            obj.alpha_ref = obj.W_alpha(:, 1);
            obj.w_ref = zeros(3,1);
        end
        
        function ref = getCommand(obj, r_atual)
            wp_alvo = obj.W_r(:, obj.idx);
            alpha_alvo = obj.W_alpha(:, obj.idx);
            
            % 1. LÓGICA DE TROCA DE WAYPOINT
            % A troca ocorre quando o drone REAL chega perto do alvo
            if norm(r_atual - wp_alvo) < obj.R_acc && obj.idx < obj.n_w
                obj.idx = obj.idx + 1;
                wp_alvo = obj.W_r(:, obj.idx);
                alpha_alvo = obj.W_alpha(:, obj.idx);
            end
            
            % =========================================================
            % 2. MALHA PD DE POSIÇÃO (O DRONE VIRTUAL / FANTASMA)
            % =========================================================
            % Vetor erro entre o Fantasma e o Waypoint
            e_p = wp_alvo - obj.r_ref;
            
            % Velocidade desejada do Fantasma (Ação Proporcional)
            v_des = obj.wn_ref * e_p; 
            if norm(v_des) > obj.v_max
                v_des = obj.v_max * (v_des / norm(v_des)); % Saturação Esférica
            end
            
            % Aceleração desejada do Fantasma (Ação Derivativa)
            Kd = obj.wn_ref * 2.0; % Amortecimento crítico geométrico
            a_des = Kd * (v_des - obj.v_ref);
            if norm(a_des) > obj.a_max
                a_des = obj.a_max * (a_des / norm(a_des)); % Saturação Esférica
            end
            
            % Integração exata da Cinemática
            obj.r_ref = obj.r_ref + obj.v_ref * obj.Ts + 0.5 * a_des * (obj.Ts^2);
            obj.v_ref = obj.v_ref + a_des * obj.Ts;
            obj.a_ref = a_des;
            
            % =========================================================
            % 3. MALHA PD DE ATITUDE (Evita o pico de -1700 Nm)
            % =========================================================
            % Erro de ângulo sempre pelo caminho mais curto (-pi a pi)
            erro_alpha = alpha_alvo - obj.alpha_ref;
            erro_alpha = atan2(sin(erro_alpha), cos(erro_alpha)); 
            
            % Ganhos da atitude virtual (geralmente mais rápidos que posição)
            Kp_att = obj.wn_ref * 3.0;
            Kd_att = obj.wn_ref * 4.0;
            
            % Taxa angular desejada
            w_des = Kp_att * erro_alpha;
            w_max = 30 * (pi/180); % Limite severo de giro: 30 graus por segundo
            if norm(w_des) > w_max
                w_des = w_max * (w_des / norm(w_des));
            end
            
            % Aceleração angular
            alpha_acc = Kd_att * (w_des - obj.w_ref);
            
            % Integração exata
            obj.alpha_ref = obj.alpha_ref + obj.w_ref * obj.Ts + 0.5 * alpha_acc * (obj.Ts^2);
            obj.w_ref = obj.w_ref + alpha_acc * obj.Ts;
            
            % =========================================================
            % 4. PACOTE DE RETORNO PARA A cControl
            % =========================================================
            ref.r_bar = obj.r_ref;
            ref.v_bar = obj.v_ref;
            ref.a_bar = obj.a_ref;
            ref.alpha_bar = obj.alpha_ref;
            ref.w_bar = obj.w_ref;
        end
    end
end