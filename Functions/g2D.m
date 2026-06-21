% g2D
% Descrição: Conversão de vetor de Gibbs para matriz de cossenos diretores

function D = g2D(g)
    D = (1+g'*g)*((1+g'*g)*eye(3)+2*(g*g')-2*skew(g));
end