clear
clc

% Calcola la funzione che descrive la mean camber line di un profilo alare
% e la sua derivata prima rispetto a x. Questa funzione calcola la linea
% media segmentando la superficie del profilo e tracciando la linea che
% congiunge i punti medi dei segmenti che compongono il profilo 

% IMPORTANTE: 
% funziona solo se le coppie di coordinate che descrivono il profilo sono
% divise fra dorso e ventre, prima dorso (y>0) poi ventre (y<0)

%% CALCOLO DELLA MEAN CAMBER LINE

% Dati del profilo alare
if ismac
    A = importdata('Theodorsen\k2.txt');
else
    A = importdata('Theodorsen/k2.txt');
end

A = getfield(A,"data"); % Matrice delle coppie x y
x = A(:,1); % Ascisse dei punti del profilo
y = A(:,2); % Ordinate dei punti del profilo

% Separazione di dorso e ventre
i = 1;
while A(i,2) ~= 0
    x_dorso(i) = A(i,1);
    y_dorso(i) = A(i,2);
    i = i + 1;
end

k = 0;
for k = 1:(length(A)-length(x_dorso))
    x_ventre(k) = A(k+length(x_dorso),1);
    y_ventre(k) = A(k+length(x_dorso),2);
end

% Interpolazione per ottenere la camber line
n = 1e6;
% t = linspace(0,pi,n);
% x_interp = 0.5*(cos(t) + 1);
x_interp = linspace(0,1,n); % Stabilisco le x a cui corrisponderanno le y di dorso e ventre
y_dorso_interp = spline(x_dorso,y_dorso,x_interp);
y_ventre_interp = spline(x_ventre,y_ventre,x_interp);
y_lm = (y_dorso_interp + y_ventre_interp) / 2;
% y_lm(1)= 0; 

% Calcolo piccole deviazioni della linea media per calcolare la derivata
h = 1e-6;
x_interp_der = linspace(0,1,n) + h; % Punti comuni
y_dorso_interp = spline(x_dorso,y_dorso,x_interp_der);
y_ventre_interp = spline(x_ventre,y_ventre,x_interp_der);
dy_lm = (y_dorso_interp + y_ventre_interp) / 2;

% Calcolo delle derivate
dy_dx = (dy_lm - y_lm) ./ h;

% Grafico
figure;
plot(x, y, 'b-', 'DisplayName', 'Airfoil geometry');
axis equal
hold on;
plot(x_interp, y_lm, 'r-', 'DisplayName', 'Mean camber line');
hold on
plot(x_interp,dy_dx,'DisplayName','Derivative of the mcl');
xlabel('x'); ylabel('y');
legend;
grid on;

%% CALCOLO DELL'ANGOLO DI THEODORSEN

eps = 1e-3; % la scelta di epsilon viene fatta guardando il plot della derivata prima, escludendo le regioni in cui inizia a oscillare 
l = linspace(-0.5 + eps, 0.5 - eps,n); % estremi di integrazione
I_aTh = dy_dx./sqrt(0.25 - l.^2); % integrale per l'angolo di Theodorsen
I_a0 = dy_dx.*sqrt((0.5+l)./(0.5-l)); % integrale per l'angolo di zero lift

a_Th = trapz(l,I_aTh) /  pi;
a_Th = rad2deg(a_Th) % angolo d'attacco di Theodorsen = -1.2625 ; da xfoil dovrebbe venire -1.2108 : errore del 4% 

a_0 = trapz(l,I_a0) * (2/pi);
a_0 = rad2deg(a_0) % angolo d'attacco a portanza nulla = -1.6835 ; da xfoil dovrebbe venire -1.968 : errore del 14%

%%
% eta = acos(-2.*x_interp);
% dy_eta = 
% a_0 = a_Th - (1/pi) * trapz(cos(eta)*dy_eta)