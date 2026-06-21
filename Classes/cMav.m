classdef cMav < handle
    properties
        m, g, Jt, Jt_inv, Ts
        mum, km, kf, w_min, w_max, n_r, sigma, D_rb, Js, G
        rho, Aa, c, CD0, CDa, CDq, CDde, CYb, CYp, CYr, CYda, CYdr
        CL0, CLa, CLq, CLde, Clb, Clp, Clr, Clda, Cldr, Cm0, Cma, Cmq, Cmde
        Cnb, Cnp, Cnr, Cnda, Cndr
        
        r, v, q, w, varpi
        f_aero, tau_aero, tau_tilde
    end
    
    methods
        function obj = cMav(p)
            obj.m = p.m; obj.g = p.g; obj.Jt = p.Jt; obj.Jt_inv = p.Jt_inv; obj.Ts = p.Ts;
            obj.mum = p.mum; obj.km = p.km; obj.kf = p.kf; obj.w_min = p.w_min; obj.w_max = p.w_max;
            obj.n_r = p.n_r; obj.sigma = p.sigma; obj.D_rb = p.D_rb; obj.Js = p.Js; obj.G = p.G;
            obj.rho = p.rho; obj.Aa = p.Aa; obj.c = p.c;
            obj.CD0 = p.CD0; obj.CDa = p.CDa; obj.CDq = p.CDq; obj.CDde = p.CDde;
            obj.CYb = p.CYb; obj.CYp = p.CYp; obj.CYr = p.CYr; obj.CYda = p.CYda; obj.CYdr = p.CYdr;
            obj.CL0 = p.CL0; obj.CLa = p.CLa; obj.CLq = p.CLq; obj.CLde = p.CLde;
            obj.Clb = p.Clb; obj.Clp = p.Clp; obj.Clr = p.Clr; obj.Clda = p.Clda; obj.Cldr = p.Cldr;
            obj.Cm0 = p.Cm0; obj.Cma = p.Cma; obj.Cmq = p.Cmq; obj.Cmde = p.Cmde;
            obj.Cnb = p.Cnb; obj.Cnp = p.Cnp; obj.Cnr = p.Cnr; obj.Cnda = p.Cnda; obj.Cndr = p.Cndr;
            
            obj.r = zeros(3,1); obj.v = zeros(3,1); 
            obj.q = [1; 0; 0; 0]; 
            obj.w = zeros(3,1); obj.varpi = zeros(obj.n_r, 1);
            obj.f_aero = zeros(3,1); obj.tau_aero = zeros(3,1); obj.tau_tilde = zeros(3,1);
        end
        
        function updateDisturbances(obj, eta_prev, delta)
            D_bg = q2D(obj.q);
            v_b = D_bg * obj.v;
            v_norm = norm(v_b);
            
            if v_b(1) > 2.0 
                % CORREÇÃO Z-UP: V_z negativo gera alpha positivo (vento de baixo)
                alpha = max(min(atan2(-v_b(3), v_b(1)), 20*pi/180), -20*pi/180);
                beta  = max(min(asin(v_b(2) / v_norm), 10*pi/180), -10*pi/180);
                
                % CORREÇÃO MATRIZ EIXO Y: Sem o sinal de menos no alpha
                D_eb = rotx(0)*roty(alpha)*rotz(beta); 
                
                c_2v = obj.c / (2 * v_norm);
                CD = obj.CD0 + obj.CDa*alpha + obj.CDq*c_2v*obj.w(2) + obj.CDde*delta(2);
                CY = obj.CYb*beta + obj.CYp*c_2v*obj.w(1) + obj.CYr*c_2v*obj.w(3) + obj.CYda*delta(1) + obj.CYdr*delta(3);
                CL = obj.CL0 + obj.CLa*alpha + obj.CLq*c_2v*obj.w(2) + obj.CLde*delta(2);
                Cl = obj.Clb*beta + obj.Clp*c_2v*obj.w(1) + obj.Clr*c_2v*obj.w(3) + obj.Clda*delta(1) + obj.Cldr*delta(3);
                Cm = obj.Cm0 + obj.Cma*alpha + obj.Cmq*c_2v*obj.w(2) + obj.Cmde*delta(2);
                Cn = obj.Cnb*beta + obj.Cnp*c_2v*obj.w(1) + obj.Cnr*c_2v*obj.w(3) + obj.Cnda*delta(1) + obj.Cndr*delta(3);
                q_dyn = 0.5 * obj.rho * v_norm^2;
                
                % Mantém a correção anterior do Lift positivo no Z
                obj.f_aero   = D_eb' * (q_dyn * obj.Aa * [-CD; CY; CL]);
                obj.tau_aero = q_dyn * obj.Aa * obj.c * [Cl; Cm; Cn];
            else
                obj.f_aero = zeros(3,1); obj.tau_aero = zeros(3,1);
            end
            
            varpi_dot = -(1 ./ obj.mum) .* obj.varpi + (obj.km ./ obj.mum) .* eta_prev;
            h_spin = zeros(3,1); h_spin_dot = zeros(3,1);
            for i = 1:obj.n_r
                J_rot = obj.D_rb(:,:,i) * obj.Js(:,:,i); 
                h_spin     = h_spin + J_rot * (obj.sigma(i) * obj.varpi(i) * [0;0;1]);
                h_spin_dot = h_spin_dot + J_rot * (obj.sigma(i) * varpi_dot(i) * [0;0;1]);
            end
            obj.tau_tilde = cross(h_spin, obj.w) - h_spin_dot;
        end
        
        function obj = integrate(obj, eta, ~)
            x = [obj.r; obj.v; obj.q; obj.w; obj.varpi]; 
            
            k1 = obj.Ts * fun(obj, x, eta);
            k2 = obj.Ts * fun(obj, x + k1/2, eta);
            k3 = obj.Ts * fun(obj, x + k2/2, eta);            
            k4 = obj.Ts * fun(obj, x + k3, eta); 
            
            xn = x + k1/6 + k2/3 + k3/3 + k4/6;
            
            obj.r = xn(1:3); obj.v = xn(4:6);
            obj.q = xn(7:10) / norm(xn(7:10)); 
            obj.w = xn(11:13);
            obj.varpi = min(max(xn(14:end), obj.w_min), obj.w_max);
        end
        
        function xp = fun(obj, x, eta)
            v1 = x(4:6); q1 = x(7:10); w1 = x(11:13); varpi1 = x(14:end);
            D1 = q2D(q1);
            
            varpi_dot = -(1 ./ obj.mum) .* varpi1 + (obj.km ./ obj.mum) .* eta;
            
            f_rot = obj.kf .* sign(varpi1) .* (varpi1.^2);
            nu = obj.G * f_rot;
            f_b = nu(1:3); tau_b = nu(4:6);
            
            r_dot = v1;
            v_dot = (1/obj.m) * (D1' * (f_b + obj.f_aero)) - [0; 0; obj.g];
            q_dot = 0.5 * [0, -w1'; w1, -skew(w1)] * q1;
            w_dot = obj.Jt_inv * (cross(obj.Jt * w1, w1) + tau_b + obj.tau_aero + obj.tau_tilde);
            
            xp = [r_dot; v_dot; q_dot; w_dot; varpi_dot];
        end
    end
end