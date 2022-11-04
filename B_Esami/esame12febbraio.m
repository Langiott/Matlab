%% esame 12 febbraio 2019
clear all; close all; clc;
%% esercizio 1 : Eulero, Runge-Kutta

f = @(t,y)-t*cos(y)+sin(t);
y0 = 0.8;

a = 0; b = 8*pi;

intervalli = [10 20 50 100] * 10;
punti = intervalli + 1;

subplot(2,1,1)
for i = 1 : length(punti)
    
    x = linspace(a,b,punti(i));
    h = x(2)-x(1);
    h2 = h/2;
    yrk = x;
    yeu = x;
    yeu = y0;
    yrk = y0;
    
    for j = 1 : length(x)-1
        
        k1 = f(x(j),yrk(j));
        k2 = f(x(j)+h2,yrk(j)+h2*k1);
        k3 = f(x(j)+h2,yrk(j)+h2*k2);
        k4 = f(x(j)+h,yrk(j)+h*k3);
        
        yrk(j+1) = yrk(j) + h*(k1+2*k2+2*k3+k4)/6;
        yeu(j+1) = yeu(j) + h*f(x(j),yeu(j));
    end
    plot(x,yrk,x,yeu)
    hold on
end
legend('RK 100', 'Eu 100','RK 200', 'Eu 200','RK 500', 'Eu 500','RK 1000', 'Eu 1000')

%% esercizio 2 : lunghezza della curva

f1 = @(x)cos(x);    % derivata della funzione ricavata dalla retta parametrica
clear f1;
fint = @(x)sqrt(1+(cos(x))^2); % funzione integranda

a = 0; b = pi;
prec = 1e-6;

n = 2;

int_value = 0; int_temp = 0;

err = inf;

while err > prec
    
    n = n+2;
    x = a;
    h = (b-a)/n;
    
    sp = 0;
    sd = 0;
    fab = fint(a)+fint(b);
    
    for i = 1 : (n-2)/2
        x = x+h;
        sd = sd + fint(x);
        x = x+h;
        sp = sp + fint(x);
    end
    x = x+h;
    sd = sd + fint(x);
    
    int_value = h*(fab+4*sd+2*sp)/3;
    
    err = abs(int_value-int_temp);
    int_temp = int_value;
       
end

fprintf('La lunghezza della curva e'' : %f\r', int_value)
fprintf('errore all''ultima iterazione : %e\r',err);
fprintf('precisione richiesta : %e\r', prec);
fprintf('Numero di punti utilizzati : %d\r\n', n+1);

%% esercizio 2 compito 2 : equazione di Duffing

fteta = @(t,y,dy,eps)-y-eps*y^3;
fy    = @(t,y,dy)dy;
dy0 = 1;
y0 = 0;

a = 0; b = 2*pi;
prec = 1e-6;

n = 1;

err = inf;
temp = 0;

while err > prec
    n = n+2;
    t = linspace(a,b,n);
    vel = t;
    pos = t;
    vel(1) = 1;
    pos(1) = 0;
    
    h = t(2)-t(1);
    h2 = h/2;
    
    for i = 1 : length(t)-1
        
        k1v = fteta(t(i),pos(i),vel(i),1);       
        k1p = fy(t(i),pos(i),vel(i));
        
        k2v = fteta(t(i)+h2,pos(i)+h2*k1p,vel(i)+h2*k1v,1);
        k2p = fy(t(i)+h2,pos(i)+h2*k1p,vel(i)+h2*k1v);
        
        k3v = fteta(t(i)+h2,pos(i)+h2*k2p,vel(i)+h2*k2v,1);
        k3p = fy(t(i)+h2,pos(i)+h2*k2p,vel(i)+h2*k2v);
        
        k4v = fteta(t(i)+h,pos(i)+h*k3p,vel(i)+h*k3v,1);
        k4p = fy(t(i)+h,pos(i)+h*k3p,vel(i)+h*k3v);
        
        vel(i+1) = vel(i) + h*(k1v+2*k2v+2*k3v+k4v)/6;
        pos(i+1) = pos(i) + h*(k1p+2*k2p+2*k3p+k4p)/6;
              
    end
    
    err = abs(pos(end)- temp);
    temp = pos(end);
    
end

subplot(2,1,2)
hold on
p1 = plot(t, pos,'LineWidth',1.5, 'DisplayName', 'pos : eps = 1');
p2 = plot(t,vel,'LineWidth',1.5, 'DisplayName', 'vel : eps = 1');

for j = 2 : 7
    vel(j,:) = vel(1,:);
    pos(j,:) = pos(1,:);
    vel(j,1) = 1;
    pos(j,1) = 0;
    
    h = t(2)-t(1);
    h2 = h/2;
    
    if j <=3
        eps = -.1*(j-1);
    else
        eps = .1*(j-1);
    end
    
    for i = 1 : length(t)-1
        
        k1v = fteta(t(i),pos(j,i),vel(j,i),1+eps);
        k1p = fy(t(i),pos(j,i),vel(j,i));
        k2v = fteta(t(i)+h2,pos(j,i)+h2*k1p,vel(j,i)+h2*k1v,1+eps);
        k2p = fy(t(i)+h2,pos(j,i)+h2*k1p,vel(j,i)+h2*k1v);
        k3v = fteta(t(i)+h2,pos(j,i)+h2*k2p,vel(j,i)+h2*k2v,1+eps);
        k3p = fy(t(i)+h2,pos(j,i)+h2*k2p,vel(j,i)+h2*k2v);
        k4v = fteta(t(i)+h,pos(j,i)+h*k3p,vel(j,i)+h*k3v,1+eps);
        k4p = fy(t(i)+h,pos(j,i)+h*k3p,vel(j,i)+h*k3v);
        
        vel(j,i+1) = vel(j,i) + h*(k1v+2*k2v+2*k3v+k4v)/6;
        pos(j,i+1) = pos(j,i) + h*(k1p+2*k2p+2*k3p+k4p)/6;
              
    end
    
    plot(t,pos(j,:),t,vel(j,:));
end

legend([p1 p2],'Location','southeast')







