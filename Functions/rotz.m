% rotz
% Descrição: Matriz de Cossenos Diretores (Passiva) em torno do eixo Zb

function R = rotz(ang)
    R = [ cos(ang), sin(ang), 0; 
         -sin(ang), cos(ang), 0; 
                 0,        0, 1];
end