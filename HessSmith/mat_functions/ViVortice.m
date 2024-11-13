function [U] = ViVortice(Centro, Estremo_1, Estremo_2, L2G_TransfMatrix, G2L_TransfMatrix)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Questa funzione permette di calcolare il vettore di velocità indotta da
% una linea di vortici. 
% Input:    Centro: Vettore colonna di coordinate del punto in cui va
%                   calcolata la velocità
%           Estremo_1/2: vettori colonna di coordinate degli estremi del
%                        pannello di vortici
%           L2G_TransfMatrix: Matrice 2x2 di trasformazione da coordinate
%                             locali a coordinate globali. Dovete
%                             calcolarla voi esternamente.
%           G2L_TransfMatrix: Matrice 2x2 di trasformazione da coordinate
%                             globali a coordinate locali. Dovete
%                             calcolarla voi esternamente.
% Output:   U : Vettore colonna contenente le componenti x ed y della
%               velocità indotta.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





% Trasformo da coordinate globali a coordinate locali
Centro = G2L_TransfMatrix * Centro;

Estremo_1 = G2L_TransfMatrix * Estremo_1;

Estremo_2 = G2L_TransfMatrix * Estremo_2;


%% Calcolo u e v in coordinate locali


% calcolo r1 (congiungente punto indotto - estremo 1)
r1 = Centro - Estremo_1;

% calcolo theta_1 (angolo che la congiungente r1 forma con l'asse x locale)
theta_1 = atan2(r1(2), r1(1));

% calcolo r1 (congiungente punto indotto - estremo 2)
r2 = Centro - Estremo_2;

% calcolo theta_2 (angolo che la congiungente r2 forma con l'asse x locale)
theta_2 = atan2(r2(2), r2(1));

% Fix in caso di auto-induzione
if (abs(theta_1)<10^(-8) && abs(theta_2)>3); theta_1=0; theta_2=pi; end
if (abs(theta_2)<10^(-8) && abs(theta_1)>3); theta_2=0; theta_1=-pi; end


% Calcolo le componenti della velocità
u = theta_2 - theta_1;
u = u / (2*pi);
v = (0.5/pi) * log(norm(r2)/norm(r1));



%% Converto da coordinate locali a coordinate globali

U = L2G_TransfMatrix * [u;v];



if abs(U(1)) < 10^(-12)
    U(1) = 0;
end
if abs(U(2)) < 10^(-12)
    U(2) = 0;
end

