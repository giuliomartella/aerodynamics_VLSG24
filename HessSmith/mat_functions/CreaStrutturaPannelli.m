function [Centro, Normale, Tangente, Estremo_1, Estremo_2, alpha, lunghezza, L2G_TransfMatrix, G2L_TransfMatrix] = CreaStrutturaPannelli(Corpo_Input)

NPannelli = length(Corpo_Input.x)-1;

Centro = zeros(NPannelli, 2);
Normale = zeros(NPannelli, 2);
Tangente = zeros(NPannelli, 2);
Estremo_1 = zeros(NPannelli, 2);
Estremo_2 = zeros(NPannelli, 2);
alpha = zeros(NPannelli, 1);
L2G_TransfMatrix = zeros(NPannelli, 2, 2);
G2L_TransfMatrix = zeros(NPannelli, 2, 2);



for i = 1:NPannelli
    
    Centro(i, 1) = (Corpo_Input.x(i) + Corpo_Input.x(i+1))/2;
    Centro(i, 2) = (Corpo_Input.y(i) + Corpo_Input.y(i+1))/2;
    
    Estremo_1(i, 1) = Corpo_Input.x(i);
    Estremo_1(i, 2) = Corpo_Input.y(i);
    
    Estremo_2(i, 1) = Corpo_Input.x(i+1);
    Estremo_2(i, 2) = Corpo_Input.y(i+1);
    
    dy = Estremo_2(i, 2) - Estremo_1(i, 2);
    dx = Estremo_2(i, 1) - Estremo_1(i, 1);
    
    angle = atan2(dy, dx);
    
    if (abs(angle)<10^(-12)); angle=0; end
    
    alpha(i) = angle;
    sinAngle = sin(angle);
    cosAngle = cos(angle);
    if(abs(sinAngle) < 10^(-12))
        sinAngle = 0;
    end
    if(abs(cosAngle) < 10^(-12))
        cosAngle = 0;
    end
    L2G_TransfMatrix(i, :, :) = [cosAngle ,  -sinAngle;
                                 sinAngle,  cosAngle];
                             
    G2L_TransfMatrix(i, :, :) = [cosAngle ,  sinAngle;
                                 -sinAngle,  cosAngle];
                             
    Normale(i, 1) = -sinAngle;
    Normale(i, 2) = cosAngle;
    
    Tangente(i, 1) = cosAngle;
    Tangente(i, 2) = sinAngle;
    
    lunghezza(i) = norm(Estremo_2(i, :) - Estremo_1(i, :));
    
end


