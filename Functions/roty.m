% roty
% Descrição: Matriz de Cossenos Diretores (Passiva) em torno do eixo Yb

function R = roty(ang)
    R = [cos(ang), 0, -sin(ang); 
                0, 1,         0; 
         sin(ang), 0,  cos(ang)];
end