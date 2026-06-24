% D2q
% Descrição: Conversão de matriz de cossenos diretores (DCM) para quatérnio
% Convenção: Escalar no final (q = [e; n], onde e é a parte vetorial e n é o escalar)

function q = D2q(D)
    tr = trace(D);
    n = 0.5 * sqrt(1 + tr);
    
    % Parte vetorial e montagem do quatérnio (escalar no final)
    e = 0.25 / n * [D(2,3) - D(3,2);
                    D(3,1) - D(1,3);
                    D(1,2) - D(2,1)];
                    
    q = [e; n];
end