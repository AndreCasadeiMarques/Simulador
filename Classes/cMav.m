classdef cMav < handle
    properties
        % Físico e Inércia
        m         % massa total da aeronave [kg]
        g         % aceleração da gravidade [m/s^2]
        Jt        % matriz de inércia total (corpo + rotores) [kg.m^2]
        Jt_inv    % inversa da matriz de inércia total
        
        % Simulação
        Ts        % passo de integração temporal / tempo de amostragem [s]
        
        % Sistema de Propulsão (Rotores e Motores)
        n_r       % número de rotores
        mum       % constante de tempo dos motores (\mu_m)
        km        % ganho dos motores (k_m)
        kf        % coeficiente de empuxo dos rotores (k_f)
        w_min     % velocidade angular mínima permitida nos rotores [rad/s]
        w_max     % velocidade angular máxima permitida nos rotores [rad/s]
        sigma     % sentido de rotação de cada rotor (+1 ou -1)
        D_rb      % matrizes de rotação nominais de cada rotor para o corpo
        Js        % matrizes de inércia individuais de cada rotor
        G         % matriz de alocação de controle (G)
        
        % Aerodinâmica - Propriedades Gerais
        rho       % densidade do ar (\rho) [kg/m^3]
        Aa        % área de referência aerodinâmica [m^2]
        c         % corda aerodinâmica de referência [m]
        
        % Aerodinâmica - Coeficientes de Força
        CD0; CDa; CDq; CDde;           % coeficientes de arrasto (Drag - Eixo X)
        CYb; CYp; CYr; CYda; CYdr;     % coeficientes de força lateral (Lateral - Eixo Y)
        CL0; CLa; CLq; CLde;           % coeficientes de sustentação (Lift - Eixo Z)
        
        % Aerodinâmica - Coeficientes de Momento
        Clb; Clp; Clr; Clda; Cldr;     % coeficientes de momento de rolagem (Roll)
        Cm0; Cma; Cmq; Cmde;           % coeficientes de momento de arfagem (Pitch)
        Cnb; Cnp; Cnr; Cnda; Cndr;     % coeficientes de momento de guinada (Yaw)
        
        % Variáveis de Estado
        r         % posição no referencial inercial 3x1 [m]
        v         % velocidade linear 3x1 [m/s]
        q         % atitude representada por quatérnio 4x1 [q0; q1; q2; q3]
        w         % velocidade angular da aeronave (\omega) 3x1 [rad/s]
        varpi     % velocidades angulares atuais dos rotores (\varpi) [rad/s]
        
        % Variáveis Computadas (Esforços Externos e Inerciais)
        f_aero    % força aerodinâmica resultante no referencial do corpo 3x1 [N]
        tau_aero  % torque aerodinâmico resultante no referencial do corpo 3x1 [N.m]
        tau_tilde % momento giroscópico resultante da inércia dos rotores 3x1 [N.m]
    end
    
    methods        
        %% CONSTRUTORA
        % Inicializa o objeto do MAV carregando todos os parâmetros físicos, 
        % geométricos e aerodinâmicos do arquivo de configuração (struct p), 
        % além de pré-alocar os vetores de estado e esforços.
        function obj = cMav(p)
            % Parâmetros Físicos e Geométricos
            obj.m = p.m; obj.g = p.g; obj.Jt = p.Jt; obj.Jt_inv = p.Jt_inv; obj.Ts = p.Ts;
            obj.mum = p.mum; obj.km = p.km; obj.kf = p.kf; obj.w_min = p.w_min; obj.w_max = p.w_max;
            obj.n_r = p.n_r; obj.sigma = p.sigma; obj.D_rb = p.D_rb; obj.Js = p.Js; obj.G = p.G;
            
            % Parâmetros Aerodinâmicos
            obj.rho = p.rho; obj.Aa = p.Aa; obj.c = p.c;
            obj.CD0 = p.CD0; obj.CDa = p.CDa; obj.CDq = p.CDq; obj.CDde = p.CDde;
            obj.CYb = p.CYb; obj.CYp = p.CYp; obj.CYr = p.CYr; obj.CYda = p.CYda; obj.CYdr = p.CYdr;
            obj.CL0 = p.CL0; obj.CLa = p.CLa; obj.CLq = p.CLq; obj.CLde = p.CLde;
            obj.Clb = p.Clb; obj.Clp = p.Clp; obj.Clr = p.Clr; obj.Clda = p.Clda; obj.Cldr = p.Cldr;
            obj.Cm0 = p.Cm0; obj.Cma = p.Cma; obj.Cmq = p.Cmq; obj.Cmde = p.Cmde;
            obj.Cnb = p.Cnb; obj.Cnp = p.Cnp; obj.Cnr = p.Cnr; obj.Cnda = p.Cnda; obj.Cndr = p.Cndr;
            
            % Inicialização das Variáveis de Estado
            obj.r = zeros(3,1); 
            obj.v = zeros(3,1); 
            obj.q = [1; 0; 0; 0]; 
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
            
            % --- Aerodinâmica (Avaliada apenas acima de 2.0 m/s) ---
            if v_b(1) > 2.0 
                % Correção Z-UP e extração de ângulos aerodinâmicos saturados
                alpha = max(min(atan2(-v_b(3), v_b(1)), 20*pi/180), -20*pi/180);
                beta  = max(min(asin(v_b(2) / v_norm), 10*pi/180), -10*pi/180);
                
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
                J_rot = obj.D_rb(:,:,i) * obj.Js(:,:,i); 
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
            
            % 4. Cinemática de Rotação (Derivada do Quatérnio)
            q_dot = 0.5 * [0, -w1'; w1, -skew(w1)] * q1;
            
            % 5. Dinâmica de Rotação (Equação de Euler com Efeitos Aerodinâmicos e Inerciais)
            w_dot = obj.Jt_inv * (cross(obj.Jt * w1, w1) + tau_b + obj.tau_aero + obj.tau_tilde);
            
            % Empacotamento do vetor derivada
            xp = [r_dot; v_dot; q_dot; w_dot; varpi_dot];
        end
    end
end