function posizione = verificaPosizionePunto_2(punto, planeCoeff)
    % Estrai i coefficienti del piano
    A = planeCoeff(1);
    B = planeCoeff(2);
    C = planeCoeff(3);
    D = planeCoeff(4);

    % Estrai le coordinate del punto
    x = punto(1);
    y = punto(2);
    z = punto(3);

    % Calcola Ax + By + Cz + D
    risultato = A*x + B*y + C*z + D;

    % Determina la posizione del punto rispetto al piano
    if risultato > 0
        posizione = 1;
    elseif risultato < 0
        posizione = -1;
    else
        posizione = 0;
    end
end
