clear all; clc; close all
format long

%% esercizio 1 : equazione differenziale

f = str2func(strcat('@(x,y)', '-(x*y)+sin(4*x)'));
y0 = 1;

n = [10 20 50 100]*10;

a = 0; b = 10*pi;
frk = figure();
frk.WindowState = 'maximized';
for i = 1 : length(n)
% controlliamo i vari metodi con i vari punti


    x = linspace(a,b, n(i));
    h = x(2)-x(1);
    
    len = length(x);
    
    yRK=zeros(len,1);
    yEU = zeros(len,1);
    yRK(1) = y0;
    yEU(1) = y0;
    
    for k = 1 : length(x)-1
        
        k1 = f(x(k),yRK(k));
        k2 = f(x(k)+h/2,yRK(k)+k1*h/2);
        k3 = f(x(k)+h/2,yRK(k)+k2*h/2);
        k4 = f(x(k)+h,yRK(k)+k3*h);       
        
        yRK(k+1) = yRK(k) + (k1+2*k2+2*k3+k4)*h/6;
        yEU(k+1) = yEU(k) + h*f(x(k),yEU(k)); 
    end
    
    subplot(length(n),1,i)
    plot(x, yRK,x, yEU);
    title([num2str(n(i)),' punti'])
    xlim([a b])
    legend('Runge-Kutta', 'Eulero');   
    hold on        
end

% sembra che il metodo di runge kutta diverga sempre
% non so se ho sbagliato...........

%% esercizio 2 : zeri della funzione

f2 = @(x)-x+1./(1+x.^4);
df2 = @(x)((-1-5*x^4)*(1+x^4)-(1-x-x^5)*(4*x^3))/((1+x^4)^2);

a1 = 0; b1 = 1;

xn = linspace(a1,b1,100);
y = f2(xn);
figure()
plot(xn,y)
yline(0);
hold on

% bisezione
err = inf;
prec = 1e-5;
a1 = 0.6;
counter1 = 0;

% Newton Raphson
xinit = 0.6;
err2 = inf;
counter2 = 0;

while err > prec
    
    % bisezione
    fa = f2(a1);
    fb = f2(b1);
    
    c = (a1+b1)/2;
    fc = f2(c);
    
    if(fa*fc < 0)
        b1 = c;
    elseif (fb*fc <0)
        a1 = c;
    else
        break
    end     
    
    err = abs(f2(c)-0);
    counter1 = counter1+1;
    
    % Newton Raphson
    if err2 > prec
        xinit = xinit-f2(xinit)/df2(xinit);
    
        err2=abs(f2(xinit)-0);  
        counter2 = counter2 + 1;
    end
    
    differr(counter1) = abs(err-err2);

end

zeroy = f2(c);
plot(c, zeroy, 'ro')

fprintf('zero con bisezione :\r\tx = %f\r\ty = %e\r\toperazioni : %d\r\n',c,zeroy,counter1);
fprintf('zero con Newton Raphson :\r\tx = %f\r\ty = %e\r\toperazioni : %d\r\n',xinit,f2(xinit),counter2);
plot(xinit, f2(xinit),'go');

fprintf('differenze degli errori : \r\n');
for i = 1 : length(differr)
    fprintf('\tciclo : %d\tdifferenza : %e\r', i, differr(i));
    
end