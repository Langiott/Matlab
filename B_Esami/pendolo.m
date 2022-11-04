clear all; clc; close all
format long

%% esercizio 2 secondo compito : equazioni del pendolo 

% d''x/dt = -sin(x)                         dy/dt = -sin(x) : y(0) = 1;
% x(0) = 0              = (dx/dt = y) =>    dx/dt = y       : x(0) = 0; 
% dx/dt(0) = 1          
%                   trasformiamo quindi un problema il valori iniziali del
%                   II ordine in un sistema di equazioni del I ordine e
%                   risolviamo con i metodi numerici classici (R-K IV)

f = @(t,u,v)v;
g = @(t,u,v)-sin(u);

a = 0; b = pi/2;
b = 10;

n = 10;

t = linspace(a, b, n);

h = t(2)-t(1);
h2 = 0.5*h;

u(1) = 0; v(1) = 1;

prec = 1e-6;
err = inf;

temp = 0;

while 1
    
    
    for i = 1 : n-1
    
        k1 = f(t(i), u(i), v(i));
        l1 = g(t(i), u(i), v(i));
    
        k2 = f(t(i)+h2, u(i)+k1*h2, v(i)+l1*h2);
        l2 = g(t(i)+h2, u(i)+k1*h2, v(i)+l1*h2);
    
        k3 = f(t(i)+h2, u(i)+k2*h2, v(i)+l2*h2);
        l3 = g(t(i)+h2, u(i)+k2*h2, v(i)+l2*h2);
    
        k4 = f(t(i)+h, u(i)+k3*h, v(i)+l3*h);
        l4 = g(t(i)+h, u(i)+k3*h, v(i)+l3*h);
    
        u(i+1) = u(i) + (h/6)*(k1 + 2*k2 + 2*k3 + k4);
        v(i+1) = v(i) + (h/6)*(l1 + 2*l2 + 2*l3 + l4);
    
    end    
    
    if err <= prec
        break
    else 
        n = n+2;  
        err = abs(v(end)-temp);
        temp = v(end);
        t = linspace(a, b, n);

        h = t(2)-t(1);
        h2 = 0.5*h;

        u(1) = 0; v(1) = 1;
    end    
       
end

fprintf('numero di intervalli : %d\r', n);
fprintf('errore               : %e\r\n', err);


subplot(2,1,1)
plot(t, u, 'r', t, v, 'b')
subplot(2,1,2)
plot3(t,v, u)

