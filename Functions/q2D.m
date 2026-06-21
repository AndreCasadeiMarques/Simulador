% q2D
% Descrição: Conversão de quatérnio para matriz de cossenos diretores

function D = q2D(q)
    % q = [q0; q1; q2; q3] (Escalar primeiro)
    eta = q(1); eps = q(2:4);
    % Conversão para DCM
    D = (eta^2 - eps'*eps)*eye(3) + 2*(eps*eps') - 2*eta*skew(eps);
end