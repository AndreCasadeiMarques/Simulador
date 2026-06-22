classdef cGuidance < handle
    properties
        % Configuração da Missão (Waypoints)
        W_r         % matriz de waypoints de posição 3xN [m]
        W_alpha     % matriz de waypoints de atitude (Euler) 3xN [rad]
        n_w         % número total de waypoints na missão
        idx         % índice do waypoint alvo atual
        
        % Limites de Desempenho e Sintonia
        R_acc       % raio de aceitação (tolerância) para troca de waypoint [m]
        v_max       % limite máximo de velocidade linear (saturação) [m/s]
        a_max       % limite máximo de aceleração linear (saturação) [m/s^2]
        wn_ref      % frequência natural da malha cinemática do drone virtual [rad/s]
        
        % Simulação
        Ts          % passo de integração temporal [s]
        
        % Estados do Drone Virtual (Referência Cinemática Contínua)
        r_ref; v_ref; a_ref;   % posição [m], velocidade [m/s] e aceleração [m/s^2] de referência
        alpha_ref; w_ref;      % atitude (Euler) [rad] e taxa angular [rad/s] de referência
    end
    
    methods    
        %% CONSTRUTORA
        % Inicializa o gerador de trajetórias (guiamento), carregando a missão e
        % posicionando o "drone virtual" (referência) exatamente no primeiro waypoint.
        function obj = cGuidance(p)
            % Parâmetros e Limites
            obj.W_r     = p.W_r; 
            obj.W_alpha = p.W_alpha;
            obj.R_acc   = p.R_acc; 
            obj.v_max   = p.v_max; 
            obj.a_max   = p.a_max;
            obj.wn_ref  = p.wn_ref; 
            obj.Ts      = p.Ts;
            
            % Controle de Missão
            obj.idx = 1; 
            obj.n_w = size(obj.W_r, 2);
            
            % Inicialização dos Estados do Drone Virtual
            obj.r_ref = obj.W_r(:, 1); 
            obj.v_ref = zeros(3,1); 
            obj.a_ref = zeros(3,1);
            
            obj.alpha_ref = obj.W_alpha(:, 1);
            obj.w_ref     = zeros(3,1);
        end
        
        
        %% GERAÇÃO DE COMANDOS E TRAJETÓRIA (DRONE VIRTUAL)
        % Executa a malha cinemática do guiamento, atualizando a posição do alvo 
        % (troca de waypoint) e integrando as equações de movimento do drone 
        % virtual para garantir uma referência de voo suave e contínua.
        function ref = getCommand(obj, r_atual)
            % Extração do Alvo Atual
            wp_alvo    = obj.W_r(:, obj.idx);
            alpha_alvo = obj.W_alpha(:, obj.idx);
            
            % 1. Lógica de Troca de Waypoint (Máquina de Estados)
            % A troca ocorre apenas quando o drone REAL entra na esfera de aceitação do alvo
            if norm(r_atual - wp_alvo) < obj.R_acc && obj.idx < obj.n_w
                obj.idx = obj.idx + 1;
                wp_alvo = obj.W_r(:, obj.idx);
                alpha_alvo = obj.W_alpha(:, obj.idx);
            end
            
            % --- Malha Cinemática de Posição ---
            
            % 2. Ação Proporcional (Velocidade Desejada) com Saturação
            e_p = wp_alvo - obj.r_ref; % Erro entre o Fantasma e o Waypoint
            
            v_des = obj.wn_ref * e_p; 
            if norm(v_des) > obj.v_max
                v_des = obj.v_max * (v_des / norm(v_des)); % Saturação Esférica
            end
            
            % 3. Ação Derivativa (Aceleração Desejada) com Saturação
            Kd = obj.wn_ref * 2.0; % Amortecimento crítico geométrico
            a_des = Kd * (v_des - obj.v_ref);
            if norm(a_des) > obj.a_max
                a_des = obj.a_max * (a_des / norm(a_des)); % Saturação Esférica
            end
            
            % 4. Integração Numérica (Translação do Drone Virtual)
            obj.r_ref = obj.r_ref + obj.v_ref * obj.Ts + 0.5 * a_des * (obj.Ts^2);
            obj.v_ref = obj.v_ref + a_des * obj.Ts;
            obj.a_ref = a_des;
            
            % --- Malha Cinemática de Atitude ---
            
            % 5. Ação Proporcional e Derivativa de Rotação
            erro_alpha = alpha_alvo - obj.alpha_ref;
            erro_alpha = atan2(sin(erro_alpha), cos(erro_alpha)); % Ajuste de caminho mais curto (-pi a pi)
            
            Kp_att = obj.wn_ref * 3.0; % Ganhos de atitude virtual (mais rápidos que posição)
            Kd_att = obj.wn_ref * 4.0;
            
            w_des = Kp_att * erro_alpha;
            w_max = 15 * (pi/180); % Limite severo de giro: 30 graus por segundo [rad/s]
            if norm(w_des) > w_max
                w_des = w_max * (w_des / norm(w_des));
            end
            
            alpha_acc = Kd_att * (w_des - obj.w_ref);
            
            % 6. Integração Numérica (Rotação do Drone Virtual)
            obj.alpha_ref = obj.alpha_ref + obj.w_ref * obj.Ts + 0.5 * alpha_acc * (obj.Ts^2);
            obj.w_ref = obj.w_ref + alpha_acc * obj.Ts;
            
            % --- Encapsulamento ---
            
            % 7. Pacote de Retorno para a Controladora Principal
            ref.r_bar     = obj.r_ref;
            ref.v_bar     = obj.v_ref;
            ref.a_bar     = obj.a_ref;
            ref.alpha_bar = obj.alpha_ref;
            ref.w_bar     = obj.w_ref;
        end
        
    end
end