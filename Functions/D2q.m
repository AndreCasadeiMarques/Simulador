% D2q
% Descrição: Conversão de matriz de cossenos diretores para quatérnio

function q = D2q(D)
    % Componente escalar
    trD = trace(D);
    eta = 0.5*sqrt(1+trD);
    % Componente vetorial
    epsilon = 0.25/eta*[D(2,3)-D(3,2);
                        D(3,1)-D(1,3);
                        D(1,2)-D(2,1)];
    % Expressão do quatérnio
    q = [eta;epsilon];
end