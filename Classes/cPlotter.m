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

            disp('>> Gráficos exibidos interativamente e PDFs salvos na pasta Outputs.');
        end
    end
end