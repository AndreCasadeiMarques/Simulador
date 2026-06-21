% tr
% Descrição: Traço de uma matriz

function trD = tr(D)
    trD = 0;
    for i = 1:3
        trD = trD + D(i,i);
    end
end