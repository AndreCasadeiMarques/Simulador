classdef cMav < handle
    properties
        % Parâmetros Físicos e Inerciais
        m         % massa total da aeronave [kg]
        g         % aceleração da gravidade [m/s^2]
        Jt        % matriz de inércia total (corpo + rotores) [kg.m^2]
        Jt_inv    % inversa da matriz de inércia total [kg.m^2]
        Ts        % passo de tempo da simulação [s]
        
        % Parâmetros dos Atuadores (Motores e Rotores)
        mum       % constante de tempo dos motores [s] (vetor n_r x 1)
        km        % ganho dos motores [rad/s] (vetor n_r x 1)
        kf        % coeficiente de empuxo dos rotores (vetor n_r x 1)
        w_min     % velocidade angular mínima do rotor (vetor n_r x 1)
        w_max     % velocidade angular máxima do rotor (vetor n_r x 1)
        n_r       % número total de rotores
        sigma     % sentido de rotação dos rotores (+1 ou -1)
        D_rb      % matrizes de rotação local dos rotores em relação ao corpo (3x3xN)
        Js        % tensores de inércia locais dos rotores (3x3xN)
        G         % matriz de alocação geométrica global
        
        % Parâmetros Aerodinâmicos
        rho; Aa; c;
        CD0; CDa; CDq; CDde;
        CYb; CYp; CYr; CYda; CYdr;
        CL0; CLa; CLq; CLde;
        Clb; Clp; Clr; Clda; Cldr;
        Cm0; Cma; Cmq; Cmde;
        Cnb; Cnp; Cnr; Cnda; Cndr;
        
        % Limites e Parâmetros de Diagnóstico/Física
        v_aero_min; v_diag_min; alpha_max; beta_max; ground_friction; gamma_v;
        
        % Variáveis de Estado
        r         % posição no referencial inercial 3x1 [m]
        v         % velocidade linear 3x1 [m/s]
        q         % atitude representada por quatérnio 4x1 [e; n] (escalar no final)
        w         % velocidade angular da aeronave (omega) 3x1 [rad/s]
        varpi     % velocidades angulares atuais dos rotores (varpi) [rad/s]
        
        % Variáveis Computadas (Esforços Externos e Inerciais)
        f_aero    % força aerodinâmica resultante no referencial do corpo 3x1 [N]
        tau_aero  % torque aerodinâmico resultante no referencial do corpo 3x1 [N.m]
        tau_tilde % momento giroscópico resultante da inércia dos rotores 3x1 [N.m]
    end
    
    methods        
        %% CONSTRUTORA
        % Inicializa o objeto do MAV carregando todos os parâmetros físicos,
        % geométricos e aerodinâmicos do arquivo de configuração (struct sMav),
        % além de pré-alocar os vetores de estado e esforços.
        function obj = cMav(sMav)
            % Parâmetros Físicos e Geométricos
            obj.m = sMav.m; obj.g = sMav.g; obj.Jt = sMav.Jt; obj.Jt_inv = sMav.Jt_inv; obj.Ts = sMav.Ts;
            obj.mum = sMav.mum; obj.km = sMav.km; obj.kf = sMav.kf; obj.w_min = sMav.w_min; obj.w_max = sMav.w_max;
            obj.n_r = sMav.n_r; obj.sigma = sMav.sigma; obj.D_rb = sMav.D_rb; obj.Js = sMav.Js; obj.G = sMav.G;
            
            % Parâmetros Aerodinâmicos
            obj.rho = sMav.rho; obj.Aa = sMav.Aa; obj.c = sMav.c;
            obj.CD0 = sMav.CD0; obj.CDa = sMav.CDa; obj.CDq = sMav.CDq; obj.CDde = sMav.CDde;
            obj.CYb = sMav.CYb; obj.CYp = sMav.CYp; obj.CYr = sMav.CYr; obj.CYda = sMav.CYda; obj.CYdr = sMav.CYdr;
            obj.CL0 = sMav.CL0; obj.CLa = sMav.CLa; obj.CLq = sMav.CLq; obj.CLde = sMav.CLde;
            obj.Clb = sMav.Clb; obj.Clp = sMav.Clp; obj.Clr = sMav.Clr; obj.Clda = sMav.Clda; obj.Cldr = sMav.Cldr;
            obj.Cm0 = sMav.Cm0; obj.Cma = sMav.Cma; obj.Cmq = sMav.Cmq; obj.Cmde = sMav.Cmde;
            obj.Cnb = sMav.Cnb; obj.Cnp = sMav.Cnp; obj.Cnr = sMav.Cnr; obj.Cnda = sMav.Cnda; obj.Cndr = sMav.Cndr;
            
            % Limites e Parâmetros Especiais
            obj.v_aero_min = sMav.v_aero_min; obj.v_diag_min = sMav.v_diag_min;
            obj.alpha_max = sMav.alpha_max; obj.beta_max = sMav.beta_max;
            obj.ground_friction = sMav.ground_friction; obj.gamma_v = sMav.gamma_v;
            
            % Inicialização das Variáveis de Estado
            obj.r = zeros(3,1); 
            obj.v = zeros(3,1); 
            obj.q = [0; 0; 0; 1];  % Inicialização com quatérnio identidade (escalar no final)
            obj.w = zeros(3,1); 
            obj.varpi = zeros(obj.n_r, 1);
            
            % Inicialização dos Esforços Externos
            obj.f_aero = zeros(3,1); 
            obj.tau_aero = zeros(3,1); 
            obj.tau_tilde = zeros(3,1);
        end
        
        
        %% ESFORÇOS AERODINÂMICOS E DISTÚRBIOS INERCIAIS
        % Atualiza os ângulos de ataque e derrapagem, extrai as forças e
        % momentos aerodinâmicos baseados na pressão dinâmica, e calcula o
        % torque giroscópico reativo gerado pela inércia de rotação dos rotores.
        function updateDisturbances(obj, eta_prev, delta)
            D_bg = q2D(obj.q);
            v_b = D_bg * obj.v;
            v_norm = norm(v_b);
            
            % --- Aerodinâmica ---
            if v_b(1) > obj.v_aero_min 
                % Correção Z-UP e extração de ângulos aerodinâmicos saturados
                alpha = max(min(atan2(-v_b(3), v_b(1)), obj.alpha_max), -obj.alpha_max);
                beta  = max(min(asin(v_b(2) / v_norm), obj.beta_max), -obj.beta_max);
                
                % Matriz de rotação Vento -> Corpo
                D_eb = rotx(0)*roty(alpha)*rotz(beta); 
                
                % Cálculo dos coeficientes aerodinâmicos dimensionais
                c_2v = obj.c / (2 * v_norm);
                CD = obj.CD0 + obj.CDa*alpha + obj.CDq*c_2v*obj.w(2) + obj.CDde*delta(2);
                CY = obj.CYb*beta + obj.CYp*c_2v*obj.w(1) + obj.CYr*c_2v*obj.w(3) + obj.CYda*delta(1) + obj.CYdr*delta(3);
                CL = obj.CL0 + obj.CLa*alpha + obj.CLq*c_2v*obj.w(2) + obj.CLde*delta(2);
                Cl = obj.Clb*beta + obj.Clp*c_2v*obj.w(1) + obj.Clr*c_2v*obj.w(3) + obj.Clda*delta(1) + obj.Cldr*delta(3);
                Cm = obj.Cm0 + obj.Cma*alpha + obj.Cmq*c_2v*obj.w(2) + obj.Cmde*delta(2);
                Cn = obj.Cnb*beta + obj.Cnp*c_2v*obj.w(1) + obj.Cnr*c_2v*obj.w(3) + obj.Cnda*delta(1) + obj.Cndr*delta(3);
                
                q_dyn = 0.5 * obj.rho * v_norm^2;
                
                % Esforços aerodinâmicos resultantes
                obj.f_aero   = D_eb' * (q_dyn * obj.Aa * [-CD; CY; CL]);
                obj.tau_aero = q_dyn * obj.Aa * obj.c * [Cl; Cm; Cn];
            else
                obj.f_aero = zeros(3,1); 
                obj.tau_aero = zeros(3,1);
            end
            
            % --- Efeito Giroscópico dos Rotores ---
            varpi_dot = -(1 ./ obj.mum) .* obj.varpi + (obj.km ./ obj.mum) .* eta_prev;
            h_spin = zeros(3,1); 
            h_spin_dot = zeros(3,1);
            
            for i = 1:obj.n_r
                J_rot = obj.D_rb(:,:,i) * obj.Js(:,:,i) * obj.D_rb(:,:,i)'; 
                h_spin     = h_spin + J_rot * (obj.sigma(i) * obj.varpi(i) * [0;0;1]);
                h_spin_dot = h_spin_dot + J_rot * (obj.sigma(i) * varpi_dot(i) * [0;0;1]);
            end
            
            obj.tau_tilde = cross(h_spin, obj.w) - h_spin_dot;
        end
        
        
        %% INTEGRAÇÃO NUMÉRICA (Runge-Kutta 4ª Ordem)
        % Propaga os estados da aeronave um passo no tempo (Ts) utilizando a
        % função de derivadas contínuas (fun).
        function obj = integrate(obj, eta, ~)
            % Montagem do vetor de estados atual
            x = [obj.r; obj.v; obj.q; obj.w; obj.varpi]; 
            
            % Passos do RK4
            k1 = obj.Ts * fun(obj, x, eta);
            k2 = obj.Ts * fun(obj, x + k1/2, eta);
            k3 = obj.Ts * fun(obj, x + k2/2, eta);            
            k4 = obj.Ts * fun(obj, x + k3, eta); 
            
            % Atualização final do estado
            xn = x + k1/6 + k2/3 + k3/3 + k4/6;
            
            % Desempacotamento e atualização das propriedades do objeto
            obj.r = xn(1:3); 
            obj.v = xn(4:6);
            
            % Restrição de Chão (Ground Constraint: z >= 0)
            if obj.r(3) < 0
                obj.r(3) = 0;
                if obj.v(3) < 0
                    obj.v(3) = 0; % Zera velocidade de descida se bateu no chão
                end
                
                % Atrito de Solo (Evitar deslizamento infinito no gelo virtual)
                obj.v(1:2) = obj.v(1:2) * obj.ground_friction;
            end
            
            obj.q = xn(7:10) / norm(xn(7:10)); % Normalização de segurança do quatérnio
            obj.w = xn(11:13);
            
            % Saturação dinâmica da velocidade dos rotores
            obj.varpi = min(max(xn(14:end), obj.w_min), obj.w_max);
        end
        
        
        %% EQUAÇÕES DE MOVIMENTO (Cinemática e Dinâmica)
        % Calcula as derivadas de todos os estados contínuos [r_dot, v_dot,
        % q_dot, w_dot, varpi_dot] num instante genérico de avaliação do integrador.
        function xp = fun(obj, x, eta)
            % Extração dos estados instantâneos
            v1 = x(4:6); 
            q1 = x(7:10); 
            w1 = x(11:13); 
            varpi1 = x(14:end);
            
            D1 = q2D(q1); % Matriz de atitude instantânea
            
            % 1. Dinâmica dos Motores (1ª Ordem)
            varpi_dot = -(1 ./ obj.mum) .* varpi1 + (obj.km ./ obj.mum) .* eta;
            
            % 2. Forças e Torques Gerados pela Propulsão
            f_rot = obj.kf .* sign(varpi1) .* (varpi1.^2);
            nu = obj.G * f_rot; % Alocação do esforço total
            f_b = nu(1:3); 
            tau_b = nu(4:6);
            
            % 3. Cinemática e Dinâmica de Translação
            r_dot = v1;
            v_dot = (1/obj.m) * (D1' * (f_b + obj.f_aero)) - [0; 0; obj.g];
            
            % 4. Cinemática de Rotação (Derivada do Quatérnio Escalar no Final)
            q_dot = 0.5 * [-skew(w1), w1; -w1', 0] * q1;
            
            % 5. Dinâmica de Rotação (Equação de Euler com Efeitos Aerodinâmicos e Inerciais)
            w_dot = obj.Jt_inv * (cross(obj.Jt * w1, w1) + tau_b + obj.tau_aero + obj.tau_tilde);
            
            % Empacotamento do vetor derivada
            xp = [r_dot; v_dot; q_dot; w_dot; varpi_dot];
        end
        
        
        %% HOVER INPUT
        % Calcula analiticamente a entrada de comando para a condição de hover.
        function eta_hover = getHoverInput(obj)
            f_hover   = (obj.m * obj.g) / (8 * cos(obj.gamma_v)); % Empuxo por rotor [N]
            w_hover   = sqrt(f_hover / obj.kf(1));             % Rotação de hover [rad/s]
            eta_hover = w_hover / obj.km(1);                   % Comando do motor
        end
        
        
        %% DIAGNOSTICS (Passive Data Logging)
        % Realiza todos os cálculos trigonométricos e físicos secundários 
        % para o logging de dados, mantendo a main.m livre de cálculos matemáticos.
        function diag = getDiagnostics(obj)
            D_bg = q2D(obj.q);
            v_body = D_bg * obj.v; 
            V_norm = norm(v_body);
            
            if V_norm > obj.v_diag_min
                diag.alpha = atan2(v_body(3), v_body(1));
                diag.beta  = asin(max(min(v_body(2)/V_norm, 1), -1));
            else
                diag.alpha = 0;
                diag.beta  = 0;
            end
            
            diag.pdin = 0.5 * obj.rho * V_norm^2;
            
            % Auditoria Física (Esforços de fato realizados)
            f_rotors_real = obj.kf .* (obj.varpi.^2); 
            wrench_real = obj.G * f_rotors_real; 
            diag.f_real   = wrench_real(1:3); 
            diag.tau_real = wrench_real(4:6); 
        end
    end
end