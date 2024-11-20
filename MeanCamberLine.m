clear
clc

% Calcola la funzione che descrive la mean camber line di un profilo alare
% e la sua derivata prima rispetto a x. Questa funzione calcola la linea
% media segmentando la superficie del profilo e tracciando la linea che
% congiunge i punti medi dei segmenti che compongono il profilo 

% IMPORTANTE: 
% funziona solo se le coppie di coordinate che descrivono il profilo sono
% divise fra dorso e ventre, prima dorso poi ventre

%% CALCOLO DELLA MEAN CAMBER LINE

% Dati del profilo alare
A = importdata('k2.txt');
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
n = 364;
% t = linspace(0,pi,n);
% x_interp = 0.5*(cos(t) + 1);
x_interp = linspace(0,1,n); % Stabilisco le x a cui corrisponderanno le y di dorso e ventre
y_dorso_interp = spline(x_dorso,y_dorso,x_interp);
y_ventre_interp = spline(x_ventre,y_ventre,x_interp);
y_lm = (y_dorso_interp + y_ventre_interp) / 2;
% y_lm(1)= 0; 

% Calcolo piccole deviazioni della linea media per calcolare la derivata
h = 1e-5;
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
% plot(x_interp,dy_dx,'DisplayName','Derivative of the mcl');
xlabel('x'); ylabel('y');
legend;
grid on;

%% CALCOLO DELL'ANGOLO DI THEODORSEN

l = linspace(-0.5 + n^-1, 0.5 - n^-1,n); % estremi di integrazione
I = dy_dx./sqrt(0.25 - l.^2); % integrale per l'angolo di Theodorsen
P = dy_dx./sqrt((0.5+l)./(0.5-l)); % integrale per l'angolo di zero lift

a_Th = trapz(l,I) /  pi;
a_Th = rad2deg(a_Th) % angolo d'attacco di Theodorsen 

a_0 = trapz(l,P) * (2/pi);
a_0 = rad2deg(a_0) % angolo d'attacco a portanza nulla 