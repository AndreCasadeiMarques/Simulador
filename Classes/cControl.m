classdef cControl < handle
    properties
        m, g, Jt
        K1_pos, K2_pos, K1_att, K2_att
        Theta_pinv
        f_min, f_max, kf, km
    end
    
    methods
        function obj = cControl(p)
            obj.m = p.m; obj.g = p.g; obj.Jt = p.Jt;
            obj.K1_pos = p.K1_pos; obj.K2_pos = p.K2_pos;
            obj.K1_att = p.K1_att; obj.K2_att = p.K2_att;
            obj.f_min = p.f_min; obj.f_max = p.f_max; 
            obj.kf = p.kf; obj.km = p.km;
            
            % Alocador R^5 (Ignora a força em Y, que é resolvida por Rolagem)
            %Theta = p.G([1, 3, 4, 5, 6], :); 
            obj.Theta_pinv = pinv(p.G); 
        end
        
        function [eta, log_ctrl] = compute(obj, ref, mav)
            D_bg = q2D(mav.q);
            
            % =============================================================
            % 1. Malha de Posição
            % =============================================================
            e_r = mav.r - ref.r_bar;
            e_v = mav.v - ref.v_bar;
            a_com = ref.a_bar - obj.K1_pos * e_r - obj.K2_pos * e_v;
            f_I_cmd = obj.m * (a_com + [0; 0; obj.g]) - D_bg' * mav.f_aero;
            
            % =============================================================
            % 2. Mapeamento CCA (Lógica Híbrida de Transição)
            % =============================================================
            R_z = rotz(ref.alpha_bar(3)); 
            f_psi = R_z * f_I_cmd;  
            
            Tx_req = f_psi(1); 
            Ty_req = f_psi(2);
            Tz_req = f_psi(3); 
            
            % ESTRATÉGIA HÍBRIDA LONGITUDINAL
            if Tx_req < 0
                % FREAR: O drone empina (Nose UP) e corta rotores horizontais.
                % CORREÇÃO: Restauração do sinal de MENOS (-Tx_req)
                theta_des = asin(max(min(Tx_req / max(Tz_req, 0.1), 1), -1));
                Tx_alloc  = 0; 
            else
                % ACELERAR: O drone voa nivelado e usa rotores horizontais.
                theta_des = 0; 
                Tx_alloc  = Tx_req; 
            end
            
            % Esforço Lateral sempre resolvido por Rolagem (Roll LEFT)
            phi_des = asin(max(min(-Ty_req / max(Tz_req, 0.1), 1), -1));
            
            % Matriz de Comando Passiva (Solo -> Corpo)
            D_cmd = rotx(phi_des) * roty(theta_des) * rotz(ref.alpha_bar(3));
            
            % =============================================================
            % 3. Malha de Atitude (Vetor de Gibbs)
            % =============================================================
            D_tilde = D_bg * D_cmd';
            den = max(1 + trace(D_tilde), 1e-4);
            e_alpha = (1 / den) * [D_tilde(2,3) - D_tilde(3,2);
                                   D_tilde(3,1) - D_tilde(1,3);
                                   D_tilde(1,2) - D_tilde(2,1)];
            
            alpha_com = - obj.K1_att * e_alpha - obj.K2_att * (mav.w - ref.w_bar);
            tau_b_c = obj.Jt * alpha_com + skew(mav.w) * (obj.Jt * mav.w) - mav.tau_aero;
            
            % =============================================================
            % 4. Alocação e Saturação
            % =============================================================
            u = [Tx_alloc; 0; Tz_req; tau_b_c]; 
            f_irrestrito = obj.Theta_pinv * u;
            f_star = min(max(f_irrestrito, obj.f_min), obj.f_max);
            eta = sqrt(f_star ./ (obj.kf .* (obj.km.^2)));
            
            log_ctrl.D_cmd = D_cmd;                   % Matriz de atitude comandada
            log_ctrl.f_cmd = [Tx_alloc; 0; Tz_req];   % Força comandada no corpo (3x1)
            log_ctrl.tau_cmd = tau_b_c;               % Torque comandado (3x1)
            log_ctrl.f_star = f_star;                 % Empuxo individual exigido (10x1)
        end
    end
end