classdef cPlotter < handle
    methods
        function obj = cPlotter()
        end
        
        function plotAll(obj, t, hist)
            disp('>> Gerando e exibindo gráficos interativos...');
            
            % Garante que a pasta existe antes de salvar os PDFs
            if ~exist('Outputs', 'dir')
                mkdir('Outputs');
            end
            
            % Conversão do Quatérnio armazenado para Ângulos de Euler
            N = length(t); Euler = zeros(3, N);
            for k = 1:N
                D_bg = q2D(hist.q(:,k));
                Euler(1,k) = atan2(D_bg(2,3), D_bg(3,3));
                Euler(2,k) = -asin(max(min(D_bg(1,3), 1), -1));
                Euler(3,k) = atan2(D_bg(1,2), D_bg(1,1));
            end
            Euler = Euler * (180/pi); % Radianos para Graus
            
            % -----------------------------------------------------------------
            % FIGURA 1: Trajetória 3D
            % -----------------------------------------------------------------
            f1 = figure('Name', 'Trajetória 3D');
            plot3(hist.r(1,:), hist.r(2,:), hist.r(3,:), 'b', 'LineWidth', 2); hold on;
            plot3(hist.r_bar(1,:), hist.r_bar(2,:), hist.r_bar(3,:), 'r--', 'LineWidth', 1.5);
            grid on; xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]'); 
            axis equal;
            margem = 2.0;
            xlim([min(hist.r(1,:)) - 10, max(hist.r(1,:)) + 10]);
            ylim([min(hist.r(2,:)) - 10, max(hist.r(2,:)) + 10]);
            zlim([min(hist.r(3,:)) - 2, max(hist.r(3,:)) + 2]);
            legend('Atual', 'Referência', 'Location', 'best');
            title('Rastreio de Trajetória 3D');
            exportgraphics(f1, 'Outputs/1_Trajetoria.pdf', 'ContentType', 'vector');
            
            % -----------------------------------------------------------------
            % FIGURA 2: Posição (X, Y, Z) - Estilo plot1c do Professor
            % -----------------------------------------------------------------
            f2 = figure('Name', 'Posição X, Y, Z');
            subplot(3,1,1);
            plot(t, hist.r(1,:), 'b', t, hist.r_bar(1,:), 'r--', 'LineWidth', 1.5); grid on;
            ylabel('X [m]'); legend('Real', 'Ref'); title('Posição no Tempo');
            
            subplot(3,1,2);
            plot(t, hist.r(2,:), 'b', t, hist.r_bar(2,:), 'r--', 'LineWidth', 1.5); grid on;
            ylabel('Y [m]');
            
            subplot(3,1,3);
            plot(t, hist.r(3,:), 'b', t, hist.r_bar(3,:), 'r--', 'LineWidth', 1.5); grid on;
            ylabel('Z [m]'); xlabel('Tempo [s]');
            exportgraphics(f2, 'Outputs/2_Posicao.pdf', 'ContentType', 'vector');
            
            % -----------------------------------------------------------------
            % FIGURA 3: Atitude (Euler)
            % -----------------------------------------------------------------
            f3 = figure('Name', 'Atitude (Euler)');
            plot(t, Euler(1,:), 'r', t, Euler(2,:), 'g', t, Euler(3,:), 'b', 'LineWidth', 1.5);
            grid on; ylabel('Euler [deg]'); xlabel('Tempo [s]'); 
            legend('Roll (\phi)', 'Pitch (\theta)', 'Yaw (\psi)', 'Location', 'best');
            title('Evolução da Atitude');
            exportgraphics(f3, 'Outputs/3_Atitude.pdf', 'ContentType', 'vector');

            % -----------------------------------------------------------------
            % FIGURA 4: Velocidade Linear e Angular
            % -----------------------------------------------------------------
            f4 = figure('Name', 'Velocidades');
            subplot(2,1,1);
            plot(t, hist.v(1,:), 'r', t, hist.v(2,:), 'g', t, hist.v(3,:), 'b', 'LineWidth', 1.5);
            grid on; ylabel('Linear [m/s]'); legend('v_x', 'v_y', 'v_z'); title('Velocidades do Veículo');
            
            subplot(2,1,2);
            plot(t, hist.w(1,:), 'r', t, hist.w(2,:), 'g', t, hist.w(3,:), 'b', 'LineWidth', 1.5);
            grid on; ylabel('Angular [rad/s]'); xlabel('Tempo [s]'); legend('\omega_x', '\omega_y', '\omega_z');
            exportgraphics(f4, 'Outputs/4_Velocidades.pdf', 'ContentType', 'vector');

            % -----------------------------------------------------------------
            % FIGURA 5: Esforços dos Atuadores
            % -----------------------------------------------------------------
            f5 = figure('Name', 'Dinâmica dos Rotores');
            subplot(2,1,1);
            plot(t, hist.varpi', 'LineWidth', 1.2);
            grid on; ylabel('\varpi [rad/s]'); title('Velocidade Angular dos Rotores');
            
            subplot(2,1,2);
            plot(t, hist.eta', 'LineWidth', 1.2);
            grid on; ylabel('\eta (0 a 1)'); xlabel('Tempo [s]'); title('Comando Normalizado (\eta)');
            exportgraphics(f5, 'Outputs/5_Rotores.pdf', 'ContentType', 'vector');

            %% Gráfico 6: Ângulos Aerodinâmicos (Ataque e Derrapagem)
            fig6 = figure('Name', 'Angulos Aerodinamicos', 'NumberTitle', 'off');
            subplot(2,1,1);
            plot(t, hist.alpha * (180/pi), 'b', 'LineWidth', 1.5);
            ylabel('\alpha (deg)');
            title('Ângulo de Ataque (\alpha) e Derrapagem (\beta)');
            grid on;

            subplot(2,1,2);
            plot(t, hist.beta * (180/pi), 'r', 'LineWidth', 1.5);
            ylabel('\beta (deg)');
            xlabel('Tempo (s)');
            grid on;
            saveas(fig6, 'Outputs/6_Angulos_Aerodinamicos.pdf');

            %% Gráfico 7: Forças Propulsivas (Comandadas vs Realizadas)
            fig7 = figure('Name', 'Forcas', 'NumberTitle', 'off');
            labels_f = {'F_x', 'F_y', 'F_z'};
            for i = 1:3
                subplot(3,1,i);
                plot(t, hist.f_cmd(i,:), 'k--', 'LineWidth', 1.5); hold on;
                plot(t, hist.f_real(i,:), 'b', 'LineWidth', 1.2);
                ylabel(sprintf('%s (N)', labels_f{i}));
                grid on;
                if i == 1
                    title('Forças Propulsivas no Corpo: Comandadas vs Realizadas');
                    legend('Comandada (Controlador)', 'Realizada (Motores)', 'Location', 'best');
                end
            end
            xlabel('Tempo (s)');
            saveas(fig7, 'Outputs/7_Forcas.pdf');

            %% Gráfico 8: Torques Propulsivos (Comandados vs Realizados)
            fig8 = figure('Name', 'Torques', 'NumberTitle', 'off');
            labels_tau = {'\tau_x', '\tau_y', '\tau_z'};
            for i = 1:3
                subplot(3,1,i);
                plot(t, hist.tau_cmd(i,:), 'k--', 'LineWidth', 1.5); hold on;
                plot(t, hist.tau_real(i,:), 'r', 'LineWidth', 1.2);
                ylabel(sprintf('%s (N.m)', labels_tau{i}));
                grid on;
                if i == 1
                    title('Torques Propulsivos no Corpo: Comandados vs Realizados');
                    legend('Comandado (Controlador)', 'Realizado (Motores)', 'Location', 'best');
                end
            end
            xlabel('Tempo (s)');
            saveas(fig8, 'Outputs/8_Torques.pdf');

            %% Gráfico 9: Pressão Dinâmica
            fig9 = figure('Name', 'Pressao Dinamica', 'NumberTitle', 'off');
            plot(t, hist.pdin, 'k', 'LineWidth', 1.5);
            ylabel('q (Pa)');
            xlabel('Tempo (s)');
            title('Pressão Dinâmica Durante o Voo');
            grid on;
            saveas(fig9, 'Outputs/9_Pressao_Dinamica.pdf');

            disp('>> Gráficos exibidos interativamente e PDFs salvos na pasta Outputs.');
        end
    end
end