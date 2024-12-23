
clc
close all
clear 

addpath("mat_functions\")
%% Input

TestCase = 0;

NomeProfilo = 'k2';
Chord = 1;
NPannelli = 177;

LE_X_Position = 0;
LE_Y_Position = 0;

%% Creazione profilo

% numero profilo:
% [x,y]=createProfile(NomeProfilo,NPannelli,Chord);
% Corpo.x = x ;
% Corpo.y = y ;

inputPts = importXfoilProfile(strcat(NomeProfilo, '.txt'));
inputPts = table2array(inputPts);
Corpo.x = inputPts(:,1);
Corpo.y = inputPts(:,2);

% Prima flippa i vettori
x = flipud(Corpo.x);
y = flipud(Corpo.y);

% %Rielaborazione punti - numero arbitrario di pannelli
% percentPannelliBa = 0.15; % percentuale pannelli dedicati al bordo d'attacco
% baChord = 0.02;     % inizio bordo d'attacco rispetto alla corda
% 
% %punti bordo d'attacco
% NPannelliBa = ceil(NPannelli * percentPannelliBa);
% 
% NPuntiBaD = ceil((NPannelliBa-1)/2);
% nBaD = ceil(NPannelliBa/2);
% xBaD = linspace(0,pi/2, NPuntiBaD+1);
% xBaD = xBaD(2:end);
% xBaD = fliplr(baChord * (1 - sin(xBaD)));
% 
% NPuntiBaV = floor((NPannelliBa-1)/2);
% nBaV = floor(NPannelliBa/2);
% xBaV = linspace(0,pi/2, NPuntiBaV+2);
% xBaV = xBaV(2:end-1);
% xBaV = baChord .* (1 - sin(xBaV));
% 
% 
% %punti rimanenti
% ba = find(~x);
% ba = ba(1);
% NPuntiDorso = ceil((NPannelli-NPannelliBa+2)/2);
% NPuntiVentre = floor((NPannelli-NPannelliBa+2)/2);
% xV = linspace(x(1), baChord, NPuntiVentre); 
% xD = linspace(baChord, x(end), NPuntiDorso);
% 
% 
% %riordinamento punti
% xV = [xV xBaV];
% xD = [xBaD xD];
% 
% yV = interp1([x(1:ba-1); 0], [y(1:ba-1); 0], xV);
% yD = interp1(x(ba:end), y(ba:end), xD);
% 
% x = [xV xD]';
% y = [yV yD]';


Corpo.x = x ;
Corpo.y = y ;
Corpo.x = x.*Chord;
Corpo.y = y.*Chord;

figure;
plot(x, y, 'o-')
axis equal

%% Creazione di una struttura di pannelli

[Centro, Normale, Tangente, Estremo_1, Estremo_2, ~, lunghezza, L2G_TransfMatrix, G2L_TransfMatrix] = CreaStrutturaPannelli(Corpo);
        
%% Inizializzazione matrici e vettori

% Ora che ho i pannelli, posso inizializzare la matrice ed i vettori

NCols = sum(NPannelli) + 1;
NRows = NCols;
matriceA = zeros(NRows, NCols);
TermineNoto = zeros(NRows, 1);

%% Creazione della matrice quadrata As


for i = 1:NPannelli
    index_i = i; % riga

    Centro_qui = Centro(i, :)';
    Normale_qui = Normale(i, :)';

    indexStart_colonna = 0;

        for j = 1:NPannelli
            index_j = indexStart_colonna + j;  % Colonna

            Estremo_1_qui = Estremo_1(j, :)';
            Estremo_2_qui = Estremo_2(j, :)';

            L2G_TransfMatrix_qui = squeeze(L2G_TransfMatrix(j, :, :));
            G2L_TransfMatrix_qui = squeeze(G2L_TransfMatrix(j, :, :));

            matriceA(index_i, index_j) = dot(ViSorgente(Centro_qui, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui), Normale_qui);

            matriceA(index_i, sum(NPannelli)+1) = matriceA(index_i, sum(NPannelli)+1) + dot(ViVortice(Centro_qui, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui), Normale_qui);


        end

end


%% Creazione delle componenti dei vettori a_v, c_s e c_v


Centro_Start = Centro(1, :)';
Tangente_Start = Tangente(1, :)';

Centro_End = Centro(end, :)';
Tangente_End = Tangente(end, :)';


b = 0;
for j = 1:NPannelli(1)

    index_j = j;

    Estremo_1_qui = Estremo_1(j, :)';
    Estremo_2_qui = Estremo_2(j, :)';
    L2G_TransfMatrix_qui = squeeze(L2G_TransfMatrix(j, :, :));
    G2L_TransfMatrix_qui = squeeze(G2L_TransfMatrix(j, :, :));

    a = dot(ViSorgente(Centro_Start, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui), Tangente_Start);
    b = b + dot(ViVortice(Centro_Start, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui), Tangente_Start);

    a = a + dot(ViSorgente(Centro_End, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui), Tangente_End);
    b = b + dot(ViVortice(Centro_End, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui), Tangente_End);


    matriceA(sum(NPannelli) + 1, index_j) = a;

end

matriceA(sum(NPannelli) + 1, sum(NPannelli) + 1) = b;

%% Ciclo per trovare alpha Theodorsen
alpha = -1.3:0.001:-1;
cpBa = zeros(length(alpha),1);
ba = find(~x);
ba = ba(1);

for ii=1:length(cpBa)
    al = alpha(ii);
    U_inf = 1;  % Velocità all'infinito [m/s]
    U_inf_x = U_inf * cos(deg2rad(al));
    U_inf_y = U_inf * sin(deg2rad(al));

    U_inf = [U_inf_x; U_inf_y];
    U_inf_normal = [-U_inf(2); U_inf(1)];
    U_inf_normal = U_inf_normal ./ norm(U_inf_normal);


    %% Creazione del termine noto
    for j = 1:NPannelli

        Normale_qui = Normale(j, :)';

        index = j;

        TermineNoto(index) = - dot(U_inf, Normale_qui);
    end

    Tangente_1 = Tangente(1, :)';
    Tangente_end = Tangente(end, :)';
    TermineNoto(sum(NPannelli) + 1) = - dot(U_inf, (Tangente_1 + Tangente_end));

    %% Risoluzione sistema lineare

    Soluzione = linsolve(matriceA,TermineNoto);

    %% Calcolo del cp e della velocità sui pannelli

    % compute sources induced velocity through loop
    sourceMat = zeros(NPannelli, NPannelli + 1);
    for i = 1:NPannelli
        index_i = i; % riga

        Centro_qui = Centro(i, :)';
        Tangente_qui = Tangente(i, :)';

        indexStart_colonna = 0;

            for j = 1:NPannelli
                index_j = indexStart_colonna + j;  % Colonna

                Estremo_1_qui = Estremo_1(j, :)';
                Estremo_2_qui = Estremo_2(j, :)';

                L2G_TransfMatrix_qui = squeeze(L2G_TransfMatrix(j, :, :));
                G2L_TransfMatrix_qui = squeeze(G2L_TransfMatrix(j, :, :));

                sourceMat(index_i, index_j) = dot(ViSorgente(Centro_qui, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui), Tangente_qui);
                sourceMat(index_i, sum(NPannelli)+1) = sourceMat(index_i, sum(NPannelli)+1) + dot(ViVortice(Centro_qui, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui), Tangente_qui);
            end

    end

    % compute stream velocity locally tangent
    velStreamTangent = zeros(NPannelli, 1);
    for j = 1:NPannelli

        Tangente_qui = Tangente(j, :)';
        index = j;
        velStreamTangent(index) = dot(U_inf, Tangente_qui);
    end


    % sum velocity 
    velPannelli = sourceMat * Soluzione + velStreamTangent;

    cp = 1 - abs(velPannelli) .^ 2;
    cpBa(ii) = cp(ba);
    

end

[~,ind] = min(abs(cpBa-1));
alphaTheodorsen = alpha(ind);


%% Recupero dati con alpha theodorsen
U_inf = 1;  % Velocità all'infinito [m/s]
alpha = alphaTheodorsen;   % Angolo di incidenza [°]
U_inf_x = U_inf * cos(deg2rad(alpha));
U_inf_y = U_inf * sin(deg2rad(alpha));

U_inf = [U_inf_x; U_inf_y];
U_inf_normal = [-U_inf(2); U_inf(1)];
U_inf_normal = U_inf_normal ./ norm(U_inf_normal);


NCols = sum(NPannelli) + 1;
NRows = NCols;
matriceA = zeros(NRows, NCols);
TermineNoto = zeros(NRows, 1);

%% Creazione della matrice quadrata As


for i = 1:NPannelli
    index_i = i; % riga

    Centro_qui = Centro(i, :)';
    Normale_qui = Normale(i, :)';

    indexStart_colonna = 0;

        for j = 1:NPannelli
            index_j = indexStart_colonna + j;  % Colonna

            Estremo_1_qui = Estremo_1(j, :)';
            Estremo_2_qui = Estremo_2(j, :)';

            L2G_TransfMatrix_qui = squeeze(L2G_TransfMatrix(j, :, :));
            G2L_TransfMatrix_qui = squeeze(G2L_TransfMatrix(j, :, :));

            matriceA(index_i, index_j) = dot(ViSorgente(Centro_qui, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui), Normale_qui);

            matriceA(index_i, sum(NPannelli)+1) = matriceA(index_i, sum(NPannelli)+1) + dot(ViVortice(Centro_qui, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui), Normale_qui);


        end

end


%% Creazione delle componenti dei vettori a_v, c_s e c_v


Centro_Start = Centro(1, :)';
Tangente_Start = Tangente(1, :)';

Centro_End = Centro(end, :)';
Tangente_End = Tangente(end, :)';


b = 0;
for j = 1:NPannelli(1)

    index_j = j;

    Estremo_1_qui = Estremo_1(j, :)';
    Estremo_2_qui = Estremo_2(j, :)';
    L2G_TransfMatrix_qui = squeeze(L2G_TransfMatrix(j, :, :));
    G2L_TransfMatrix_qui = squeeze(G2L_TransfMatrix(j, :, :));

    a = dot(ViSorgente(Centro_Start, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui), Tangente_Start);
    b = b + dot(ViVortice(Centro_Start, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui), Tangente_Start);

    a = a + dot(ViSorgente(Centro_End, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui), Tangente_End);
    b = b + dot(ViVortice(Centro_End, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui), Tangente_End);


    matriceA(sum(NPannelli) + 1, index_j) = a;

end

matriceA(sum(NPannelli) + 1, sum(NPannelli) + 1) = b;



%% Creazione del termine noto

for j = 1:NPannelli

    Normale_qui = Normale(j, :)';

    index = j;

    TermineNoto(index) = - dot(U_inf, Normale_qui);
end

Tangente_1 = Tangente(1, :)';
Tangente_end = Tangente(end, :)';
TermineNoto(sum(NPannelli) + 1) = - dot(U_inf, (Tangente_1 + Tangente_end));

%% Risoluzione sistema lineare

Soluzione = linsolve(matriceA,TermineNoto);

%% Calcolo del cp e della velocità sui pannelli

% compute sources induced velocity through loop
sourceMat = zeros(NPannelli, NPannelli + 1);
for i = 1:NPannelli
    index_i = i; % riga

    Centro_qui = Centro(i, :)';
    Tangente_qui = Tangente(i, :)';

    indexStart_colonna = 0;

        for j = 1:NPannelli
            index_j = indexStart_colonna + j;  % Colonna

            Estremo_1_qui = Estremo_1(j, :)';
            Estremo_2_qui = Estremo_2(j, :)';

            L2G_TransfMatrix_qui = squeeze(L2G_TransfMatrix(j, :, :));
            G2L_TransfMatrix_qui = squeeze(G2L_TransfMatrix(j, :, :));

            sourceMat(index_i, index_j) = dot(ViSorgente(Centro_qui, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui), Tangente_qui);
            sourceMat(index_i, sum(NPannelli)+1) = sourceMat(index_i, sum(NPannelli)+1) + dot(ViVortice(Centro_qui, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui), Tangente_qui);
        end

end

% compute stream velocity locally tangent
velStreamTangent = zeros(NPannelli, 1);
for j = 1:NPannelli

    Tangente_qui = Tangente(j, :)';
    index = j;
    velStreamTangent(index) = dot(U_inf, Tangente_qui);
end


% sum velocity 
velPannelli = sourceMat * Soluzione + velStreamTangent;

cp = 1 - abs(velPannelli) .^ 2;

figure()
%x = linspace(0,2*pi,50);
%y = sin(x) + randi(50,1,50);
%c = linspace(1,10,length(x));
scatter(x(1:end-1),y(1:end-1),[],-cp,'filled')
axis equal
grid on
colorbar
colormap jet
figure()
plot(x(1:end-1), -cp)
grid on


%% calcolo coeff. portanza e coeff. momento
circ=sum(velPannelli.*lunghezza(1:NPannelli)'); %circolazione
cl=zeros(1,NPannelli);
cm=zeros(1,NPannelli);      
for i=1:NPannelli
    cl(i)=cos(deg2rad(-1.984))*(-cp(i)*lunghezza(i)*Normale(i,2))-sin(deg2rad(-1.984))*(-cp(i)*lunghezza(i)*(Normale(i,1)));
end
CL=sum(cl);
CL_kj=2*circ;   %CL con formula Kutta-Joukowsky
for i=1:NPannelli
    cm(i)=cp(i)*(Centro(i,1)-0.25)*lunghezza(i)*Normale(i,2);
end
Cm=sum(cm);     %Cm rispetto al centro aerodinamico


