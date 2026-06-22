classdef cPlotter < handle
    
    methods
        
        %% CONSTRUTORA
        function obj = cPlotter()
            % A classe plotter não necessita de propriedades armazenadas,
            % atua apenas como um pacote de funções de visualização.
        end
        
        
        %% GERAÇÃO E EXPORTAÇÃO DE GRÁFICOS
        % Converte os quatérnios de atitude e plota todos os históricos de 
        % voo (trajetória, estados, atuação e aerodinâmica), salvando-os 
        % automaticamente na pasta Outputs.
        function plotAll(obj, t, hist)
            disp('>> Gerando e exibindo gráficos interativos...');
            
            % Garante que o diretório de saída existe
            if ~exist('Outputs', 'dir')
                mkdir('Outputs');
            end
            
            % --- Pré-Processamento: Quatérnio para Euler [graus] ---
            N = length(t); 
            Euler = zeros(3, N);
            
            for k = 1:N
                D_bg = q2D(hist.q(:,k));
                Euler(1,k) = atan2(D_bg(2,3), D_bg(3,3));
                Euler(2,k) = -asin(max(min(D_bg(1,3), 1), -1));
                Euler(3,k) = atan2(D_bg(1,2), D_bg(1,1));
            end
            Euler = Euler * (180/pi); 
            
            
            %% FIGURA 1: TRAJETÓRIA 3D
            f1 = figure('Name', 'Trajetória 3D', 'NumberTitle', 'off');
            
            plot3(hist.r(1,:), hist.r(2,:), hist.r(3,:), 'b', 'LineWidth', 2); hold on;
            plot3(hist.r_bar(1,:), hist.r_bar(2,:), hist.r_bar(3,:), 'r--', 'LineWidth', 1.5);
            
            grid on; axis equal;
            xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]'); 
            title('Rastreio de Trajetória 3D');
            legend('Atual', 'Referência', 'Location', 'best');
            
            % Ajuste dinâmico de margem
            xlim([min(hist.r(1,:)) - 10, max(hist.r(1,:)) + 10]);
            ylim([min(hist.r(2,:)) - 10, max(hist.r(2,:)) + 10]);
            zlim([min(hist.r(3,:)) - 2, max(hist.r(3,:)) + 2]);
            
            exportgraphics(f1, 'Outputs/1_Trajetoria.pdf', 'ContentType', 'vector');
            
            
            %% FIGURA 2: POSIÇÃO NO TEMPO (X, Y, Z)
            f2 = figure('Name', 'Posição X, Y, Z', 'NumberTitle', 'off');
            
            subplot(3,1,1);
            plot(t, hist.r(1,:), 'b', t, hist.r_bar(1,:), 'r--', 'LineWidth', 1.5); 
            grid on; ylabel('X [m]'); title('Posição no Tempo'); legend('Real', 'Ref');
            
            subplot(3,1,2);
            plot(t, hist.r(2,:), 'b', t, hist.r_bar(2,:), 'r--', 'LineWidth', 1.5); 
            grid on; ylabel('Y [m]');
            
            subplot(3,1,3);
            plot(t, hist.r(3,:), 'b', t, hist.r_bar(3,:), 'r--', 'LineWidth', 1.5); 
            grid on; ylabel('Z [m]'); xlabel('Tempo [s]');
            
            exportgraphics(f2, 'Outputs/2_Posicao.pdf', 'ContentType', 'vector');
            
            
            %% FIGURA 3: ATITUDE (ÂNGULOS DE EULER)
            f3 = figure('Name', 'Atitude (Euler)', 'NumberTitle', 'off');
            
            plot(t, Euler(1,:), 'r', t, Euler(2,:), 'g', t, Euler(3,:), 'b', 'LineWidth', 1.5);
            grid on; 
            ylabel('Euler [deg]'); xlabel('Tempo [s]'); 
            title('Evolução da Atitude');
            legend('Roll (\phi)', 'Pitch (\theta)', 'Yaw (\psi)', 'Location', 'best');
            
            exportgraphics(f3, 'Outputs/3_Atitude.pdf', 'ContentType', 'vector');
            
            
            %% FIGURA 4: VELOCIDADES LINEARES E ANGULARES
            f4 = figure('Name', 'Velocidades', 'NumberTitle', 'off');
            
            subplot(2,1,1);
            plot(t, hist.v(1,:), 'r', t, hist.v(2,:), 'g', t, hist.v(3,:), 'b', 'LineWidth', 1.5);
            grid on; ylabel('Linear [m/s]'); title('Velocidades do Veículo'); legend('v_x', 'v_y', 'v_z');
            
            subplot(2,1,2);
            plot(t, hist.w(1,:), 'r', t, hist.w(2,:), 'g', t, hist.w(3,:), 'b', 'LineWidth', 1.5);
            grid on; ylabel('Angular [rad/s]'); xlabel('Tempo [s]'); legend('\omega_x', '\omega_y', '\omega_z');
            
            exportgraphics(f4, 'Outputs/4_Velocidades.pdf', 'ContentType', 'vector');
            
            
            %% FIGURA 5: ATUAÇÃO (ROTAÇÃO E COMANDO DOS MOTORES)
            f5 = figure('Name', 'Dinâmica dos Rotores', 'NumberTitle', 'off');
            
            subplot(2,1,1);
            plot(t, hist.varpi', 'LineWidth', 1.2);
            grid on; ylabel('\varpi [rad/s]'); title('Velocidade Angular dos Rotores');
            
            subplot(2,1,2);
            plot(t, hist.eta', 'LineWidth', 1.2);
            grid on; ylabel('\eta (0 a 1)'); xlabel('Tempo [s]'); title('Comando Normalizado (\eta)');
            
            exportgraphics(f5, 'Outputs/5_Rotores.pdf', 'ContentType', 'vector');
            
            
            %% FIGURA 6: ÂNGULOS AERODINÂMICOS (ATAQUE E DERRAPAGEM)
            f6 = figure('Name', 'Angulos Aerodinamicos', 'NumberTitle', 'off');
            
            subplot(2,1,1);
            plot(t, hist.alpha * (180/pi), 'b', 'LineWidth', 1.5);
            grid on; ylabel('\alpha [deg]'); title('Ângulo de Ataque (\alpha) e Derrapagem (\beta)');
            
            subplot(2,1,2);
            plot(t, hist.beta * (180/pi), 'r', 'LineWidth', 1.5);
            grid on; ylabel('\beta [deg]'); xlabel('Tempo [s]');
            
            exportgraphics(f6, 'Outputs/6_Angulos_Aerodinamicos.pdf', 'ContentType', 'vector');
            
            
            %% FIGURA 7: AUDITORIA DE FORÇAS PROPULSIVAS (COMANDADO VS REAL)
            f7 = figure('Name', 'Forcas Propulsivas', 'NumberTitle', 'off');
            labels_f = {'F_x', 'F_y', 'F_z'};
            
            for i = 1:3
                subplot(3,1,i);
                plot(t, hist.f_cmd(i,:), 'k--', 'LineWidth', 1.5); hold on;
                plot(t, hist.f_real(i,:), 'b', 'LineWidth', 1.2);
                grid on; ylabel(sprintf('%s [N]', labels_f{i}));
                
                if i == 1
                    title('Forças Propulsivas no Corpo: Comandadas vs Realizadas');
                    legend('Comandada (Controlador)', 'Realizada (Motores)', 'Location', 'best');
                end
            end
            xlabel('Tempo [s]');
            
            exportgraphics(f7, 'Outputs/7_Forcas.pdf', 'ContentType', 'vector');
            
            
            %% FIGURA 8: AUDITORIA DE TORQUES PROPULSIVOS (COMANDADO VS REAL)
            f8 = figure('Name', 'Torques Propulsivos', 'NumberTitle', 'off');
            labels_tau = {'\tau_x', '\tau_y', '\tau_z'};
            
            for i = 1:3
                subplot(3,1,i);
                plot(t, hist.tau_cmd(i,:), 'k--', 'LineWidth', 1.5); hold on;
                plot(t, hist.tau_real(i,:), 'r', 'LineWidth', 1.2);
                grid on; ylabel(sprintf('%s [N.m]', labels_tau{i}));
                
                if i == 1
                    title('Torques Propulsivos no Corpo: Comandados vs Realizados');
                    legend('Comandado (Controlador)', 'Realizado (Motores)', 'Location', 'best');
                end
            end
            xlabel('Tempo [s]');
            
            exportgraphics(f8, 'Outputs/8_Torques.pdf', 'ContentType', 'vector');
            
            
            %% FIGURA 9: PRESSÃO DINÂMICA
            f9 = figure('Name', 'Pressao Dinamica', 'NumberTitle', 'off');
            
            plot(t, hist.pdin, 'k', 'LineWidth', 1.5);
            grid on; ylabel('q [Pa]'); xlabel('Tempo [s]');
            title('Pressão Dinâmica Durante o Voo');
            
            exportgraphics(f9, 'Outputs/9_Pressao_Dinamica.pdf', 'ContentType', 'vector');
            
            disp('>> Gráficos exibidos interativamente e PDFs salvos na pasta Outputs.');
        end
    end
end