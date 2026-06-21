% a2D
% Descrição: Conversão de ângulos de Euler para matriz de cossenos 
%            diretores (sequência 123)

function D = a2D(eul)
    % Ãngulos de Euler
    phi = eul(1);
    tht = eul(2);
    psi = eul(3);
    % Matrizes de rotação de Euler
    D3 = [cos(psi),sin(psi),0;
          -sin(psi),cos(psi),0;
          0,0,1];
    D2 = [cos(tht),0,-sin(tht);
          0,1,0;
          sin(tht),0,cos(tht)];
    D1 = [1,0,0;
          0,cos(phi),sin(phi);
          0,-sin(phi),cos(phi)];
    
    D = D3*D2*D1;
end