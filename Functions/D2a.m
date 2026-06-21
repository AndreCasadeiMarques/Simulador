% eul2D
% Descrição: Conversão matriz de cossenos diretores (sequência 123) 
%            para ângulos de Euler 

function eul = D2a(D)
    phi = -atan2(D(3,2),D(3,3));
    tht = asin(D(3,1));
    psi = -atan2(D(2,1),D(1,1));
    eul = [phi;tht;psi];
end