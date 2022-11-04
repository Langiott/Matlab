%% esame 21 aprile 2020
clc; clear all; close all;
%% Esercizio 1 : runge kutta IV stadio

f = @(x,y)(1-x)/(y^2);
fesatta = @(x)(3*x-3*x.^2/2 + 1).^(1/3);

y0 = 1;

a = 0;
b = 2;

n = 1;

while err > prec
    n = n+2;
    x = linspace(a, b, n);
    y = x;
    y(1) = y0;
    
    h = x(2)-x(1);
    h2 = h/2;
    
    for i = 1:length(x)-1
        k1 = f(x(i), y(i));
        k2 = f(x(i) + h2, y(i) + k1*h2);
        k3 = f(x(i) + h2, y(i) + k2*h2);
        k4 = f(x(i) + h, y(i) + k3*h);
        
        y(i+1) = y(i) + h*(k1 + 2*k2 + 2*k3 + k4)/6;
    end
    
    err = abs(fesatta(1)-y(x==1));
    
end

fprintf('Runge-Kutta IV stadio in x = 1 vale : %f\r', y(x==1));
fprintf('La soluzione esatta in x = 1 vale   : %f\r', fesatta(1));
fprintf('L''errore vale : %e\r', err);
fprintf('Numero di punti usati : %d\r\n', n);

yesatta = fesatta(x);
subplot(2,1,1)
plot(x, y, x, yesatta);
legend('Runge-Kutta', 'Esatta');

%% esercizio 2 : integrale con Simpson

a = -1.1*pi;
b = 0.82*pi;
b_intermedio = 0;


f0 = @(x)exp(sin(x));
f1 = @(x)-exp(-sin(x));

n = 100;

x1 = linspace(a, b_intermedio, n);
x2 = linspace(b_intermedio, b, n);


h = x1(2)-x1(1);
f1ab = f1(a) + f1(b_intermedio);
somme_pari_1 = 0; somme_dispari_1 = 0;

x = a;

while x <= b_intermedio
    
    x = x + h;
    somme_dispari_1 = somme_dispari_1 + f1(x);
    x = x + h;
    somme_pari_1 = somme_pari_1 + f1(x);
                  
end

integrale_1 = h*(f1ab + 4*somme_dispari_1 + 2*somme_pari_1)/3;

f0ab = f0(b_intermedio) + f0(b);

somme_pari_0 = 0; somme_dispari_0 = 0;

x = b_intermedio;

while x < b
    
    x = x + h;  
    somme_dispari_0 = somme_dispari_0 + f0(x);
    x = x + h;
    somme_pari_0 = somme_pari_0 + f0(x);
    
    
end

integrale_0 = h*(f0ab + 4*somme_dispari_0 + 2*somme_pari_0)/3;

integrale_totale = integrale_1 + integrale_0;

fprintf('La soluzione esatta dell''integrale : -1.0236\r')
fprintf('La soluzione con simpson e %d punti : %f\r\n',n*2, integrale_totale);

y = [f1(x1) f0(x2)];
subplot(2,1,2)
plot([x1 x2], y);














    