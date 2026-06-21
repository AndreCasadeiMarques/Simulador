% D2g
% Descrição: Conversão de matriz de cossenos diretores para vetor de Gibbs

function g = D2g(D)
    g = (1/(1+trace(D)))*[D(2,3)-D(3,2);
                          D(3,1)-D(1,3);
                          D(1,2)-D(2,1)];
end