% q2D
% Descrição: Conversão de quatérnio para matriz de cossenos diretores (DCM)
% Convenção: Escalar no final (q = [e; n], onde e é a parte vetorial e n é o escalar)

function D = q2D(q)
    e = q(1:3);
    n = q(4);
    
    % Conversão para DCM usando a formulação com escalar no final
    D = (n^2 - e'*e)*eye(3) + 2*(e*e') - 2*n*skew(e);
end