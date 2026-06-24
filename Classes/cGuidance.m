classdef cGuidance < handle
    properties
        % Configuração da Missão (Waypoints)
        W_r         % matriz de waypoints de posição 3xN [m]
        W_alpha     % matriz de waypoints de atitude (Euler) 3xN [rad]
        n_w         % número total de waypoints na missão
        idx         % índice do waypoint de origem da perna atual
        mode        % modo de voo ('armado', 'multicoptero', 'transicao')
        
        % Parâmetros (Polinômio Minimum Jerk)
        v_avg       % velocidade escalar média [m/s]
        
        % Parâmetros (Campo Vetorial - Cruzeiro)
        R_acc       % raio de wayset de cruzeiro [m]
        v_max       % velocidade máxima [m/s]
        a_max       % aceleração máxima [m/s^2]
        wn_ref      % frequência natural da malha proporcional [rad/s]
        
        % Estados do Drone Virtual (Cruzeiro)
        r_ref       
        v_ref       
        a_ref       
        alpha_ref   
        w_ref       
        
        % Estados do Polinômio (Decolagem/Pouso)
        r_0         % posição inicial [m]
        alpha_0     % atitude inicial [rad]
        t_traj      % tempo percorrido [s]
        t_j         % tempo total previsto [s]
        
        Ts          % passo de integração [s]
    end
    
    methods    
        %% CONSTRUTORA
        function obj = cGuidance(sGuidance)
            obj.W_r     = sGuidance.W_r; 
            obj.W_alpha = sGuidance.W_alpha;
            obj.v_avg   = sGuidance.v_avg; 
            obj.Ts      = sGuidance.Ts;
            obj.mode    = sGuidance.mode;
            
            obj.R_acc   = sGuidance.R_acc;
            obj.v_max   = sGuidance.v_max;
            obj.a_max   = sGuidance.a_max;
            obj.wn_ref  = sGuidance.wn_ref;
            
            obj.idx = 1; 
            obj.n_w = size(obj.W_r, 2);
            
            % Estados Iniciais
            obj.r_0       = obj.W_r(:, 1); 
            obj.alpha_0   = obj.W_alpha(:, 1);
            obj.r_ref     = obj.r_0;
            obj.v_ref     = zeros(3,1);
            obj.a_ref     = zeros(3,1);
            obj.alpha_ref = obj.alpha_0;
            obj.w_ref     = zeros(3,1);
            
            obj.t_traj  = 0.0;
            
            if obj.n_w > 1
                dist = norm(obj.W_r(:, 2) - obj.r_0);
                obj.t_j = max(dist / obj.v_avg, 0.1);
            else
                obj.t_j = 1.0; 
            end
        end
        
        %% GERAÇÃO DE COMANDOS E TRAJETÓRIA (HÍBRIDA)
        function ref = getCommand(obj, r_atual)
            
            if obj.idx == 1
                %% ---------------- FASE 1: TAKEOFF (Minimum Jerk) ----------------
                alvo_idx = 2;
                wp_alvo = obj.W_r(:, alvo_idx);
                alpha_alvo = obj.W_alpha(:, alvo_idx);
                
                obj.t_traj = obj.t_traj + obj.Ts;
                tau = min(1.0, obj.t_traj / obj.t_j);
                
                lambda = 10*tau^3 - 15*tau^4 + 6*tau^5;
                if tau < 1.0
                    lambda_dot  = (30*tau^2 - 60*tau^3 + 30*tau^4) / obj.t_j;
                    lambda_ddot = (60*tau - 180*tau^2 + 120*tau^3) / (obj.t_j^2);
                else
                    lambda_dot  = 0.0;
                    lambda_ddot = 0.0;
                end
                
                ref.r_bar = obj.r_0 + lambda * (wp_alvo - obj.r_0);
                ref.v_bar = lambda_dot * (wp_alvo - obj.r_0);
                ref.a_bar = lambda_ddot * (wp_alvo - obj.r_0);
                
                erro_alpha = alpha_alvo - obj.alpha_0;
                erro_alpha = atan2(sin(erro_alpha), cos(erro_alpha));
                ref.alpha_bar = obj.alpha_0 + lambda * erro_alpha;
                ref.w_bar     = lambda_dot * erro_alpha;
                
                % Transição para Cruzeiro (Gatilho de Tempo)
                if obj.t_traj >= obj.t_j
                    obj.idx = 2;
                    % Transfere o estado final do polinômio para o Drone Virtual
                    obj.r_ref = wp_alvo;
                    obj.v_ref = zeros(3,1);
                    obj.a_ref = zeros(3,1);
                    obj.alpha_ref = alpha_alvo;
                    obj.w_ref = zeros(3,1);
                end
                
            elseif obj.idx > 1 && obj.idx < obj.n_w - 1
                %% ---------------- FASE 2: CRUZEIRO (Campo Vetorial / Drone Virtual) ----------------
                alvo_idx = obj.idx + 1;
                wp_alvo = obj.W_r(:, alvo_idx);
                alpha_alvo = obj.W_alpha(:, alvo_idx);
                
                % Define tamanho do Wayset: "Anchoring" relaxado antes de pousar, Fly-by para os demais
                if alvo_idx == obj.n_w - 1
                    current_R_acc = 1.5; 
                else
                    current_R_acc = obj.R_acc; 
                end
                
                % Transição por entrada na esfera (Wayset)
                if norm(r_atual - wp_alvo) < current_R_acc
                    obj.idx = obj.idx + 1;
                    
                    if obj.idx == obj.n_w - 1
                        % Início da Fase de Pouso
                        obj.r_0 = wp_alvo; % Zera a velocidade matemática ancorando no waypoint
                        obj.alpha_0 = alpha_alvo;
                        obj.t_traj = 0.0;
                        dist = norm(obj.W_r(:, obj.n_w) - obj.r_0);
                        obj.t_j = max(dist / obj.v_avg, 0.1);
                        % Recalcula já no estado de pouso
                        ref = obj.getCommand(r_atual); 
                        return;
                    else
                        % Continua no cruzeiro para o próximo alvo
                        wp_alvo = obj.W_r(:, obj.idx + 1);
                        alpha_alvo = obj.W_alpha(:, obj.idx + 1);
                    end
                end
                
                % --- Cinemática do Drone Virtual (Posição) ---
                e_p = wp_alvo - obj.r_ref;
                v_des = obj.wn_ref * e_p;
                if norm(v_des) > obj.v_max
                    v_des = obj.v_max * (v_des / norm(v_des));
                end
                
                a_des = (obj.wn_ref * 2.0) * (v_des - obj.v_ref);
                if norm(a_des) > obj.a_max
                    a_des = obj.a_max * (a_des / norm(a_des));
                end
                
                obj.r_ref = obj.r_ref + obj.v_ref * obj.Ts + 0.5 * a_des * (obj.Ts^2);
                obj.v_ref = obj.v_ref + a_des * obj.Ts;
                obj.a_ref = a_des;
                
                ref.r_bar = obj.r_ref;
                ref.v_bar = obj.v_ref;
                ref.a_bar = obj.a_ref;
                
                % --- Cinemática do Drone Virtual (Atitude) ---
                % Lógica de Yaw Reativo: O drone aponta o nariz para o alvo
                vec_dir = wp_alvo(1:2) - obj.r_ref(1:2);
                if norm(vec_dir) > 0.5
                    alpha_alvo(3) = atan2(vec_dir(2), vec_dir(1));
                end
                
                erro_alpha = alpha_alvo - obj.alpha_ref;
                erro_alpha = atan2(sin(erro_alpha), cos(erro_alpha));
                
                w_des = (obj.wn_ref * 3.0) * erro_alpha;
                w_max_lim = 30 * (pi/180); % Aumentado para permitir guinadas mais ágeis
                if norm(w_des) > w_max_lim
                    w_des = w_max_lim * (w_des / norm(w_des));
                end
                
                alpha_acc = (obj.wn_ref * 4.0) * (w_des - obj.w_ref);
                
                obj.alpha_ref = obj.alpha_ref + obj.w_ref * obj.Ts + 0.5 * alpha_acc * (obj.Ts^2);
                obj.w_ref = obj.w_ref + alpha_acc * obj.Ts;
                
                ref.alpha_bar = obj.alpha_ref;
                ref.w_bar     = obj.w_ref;
                
            else
                %% ---------------- FASE 3: POUSO (Minimum Jerk) ----------------
                alvo_idx = obj.n_w;
                wp_alvo = obj.W_r(:, alvo_idx);
                alpha_alvo = obj.W_alpha(:, alvo_idx);
                
                % Lógica de "Snap" Cinemático de Chão
                if norm(r_atual - wp_alvo) < 0.2
                    obj.t_traj = obj.t_j;
                end
                
                obj.t_traj = obj.t_traj + obj.Ts;
                tau = min(1.0, obj.t_traj / obj.t_j);
                
                lambda = 10*tau^3 - 15*tau^4 + 6*tau^5;
                if tau < 1.0
                    lambda_dot  = (30*tau^2 - 60*tau^3 + 30*tau^4) / obj.t_j;
                    lambda_ddot = (60*tau - 180*tau^2 + 120*tau^3) / (obj.t_j^2);
                else
                    lambda_dot  = 0.0;
                    lambda_ddot = 0.0;
                end
                
                ref.r_bar = obj.r_0 + lambda * (wp_alvo - obj.r_0);
                ref.v_bar = lambda_dot * (wp_alvo - obj.r_0);
                ref.a_bar = lambda_ddot * (wp_alvo - obj.r_0);
                
                erro_alpha = alpha_alvo - obj.alpha_0;
                erro_alpha = atan2(sin(erro_alpha), cos(erro_alpha));
                ref.alpha_bar = obj.alpha_0 + lambda * erro_alpha;
                ref.w_bar     = lambda_dot * erro_alpha;
            end
        end
    end
end