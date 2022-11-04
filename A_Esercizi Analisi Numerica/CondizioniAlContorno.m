close all;
clear all;
clc;
n=51;% imposto il numero di passi
nm2=n-2;
xa=0;% estremo a di dove è definito il problema
xb=1;% estremo b di dove è definito il problema
ya=0.5;% sono i valori al contorno della funzione 
yb=log(2);
h=(xb-xa)/(n-1);% passo di discretizzazione 
hsq=h*h;% nella equazione compare anche hal quadrato, e me lo calcolo per sempre
x=zeros(1,n);% creo la sapazio da allocare 
p=zeros(1,n);% creo le sapazio da allocare 
q=zeros(1,n);% creo le sapazio da allocare
r=zeros(1,n);% creo le sapazio da allocare
b=zeros(n,1);
r=zeros(1,n);
x=linspace(xa,xb,n);% inizializzo definisco il vettore delle x




for i=1:n
    p(i)=fp(x(i));
    q(i)=fq(x(i));
    r(i)=fr(x(i));
end
% preparo lo spazio per costruire la matrice tridiagonale 
a=zeros(n,n);% mi fa comodo inizializzarla in questo modo in quanto mi mette gli zeri ovunque
a(1,1)=1;
a(n,n)=1;
for i=2:n-1
    a(i,i)=-2*(2+hsq*q(i));
    a(i,i-1)=2+h*p(i);
    a(i,i+1)=2-h*p(i);
end
% inizializzo il vettore dei termini noti

b(1,1)=ya;
b(n,1)=yb;
for i=2:n-1
    b(i,1)=2*hsq*r(i);
end
for i=1:n
    z(i)=zsol(x(i));% stampo la soluzione esatta 
end
% applico l'algoritmo di craut per risolvere il sistema tridiagonale
l(1,1)=a(1,1);
u(1,2)=a(1,2)/l(1,1);
y(1)=b(1)/l(1,1);
for i=2:n
    l(i,i-1)=a(i,i-1);
end
for i=2:n-1
    l(i,i-1)=a(i,i-1);
    l(i,i)=a(i,i)-l(i,i-1)*u(i-1,i);
    u(i,i+1)=a(i,i+1)/l(i,i);
    y(i)=(b(i)-l(i,i-1)*y(i-1))/l(i,i);
end
l(n,n)=a(n,n)-l(n,n-1)*u(n-1,n);
y(n)=(b(n)-l(n,n-1)*y(n-1))/l(n,n);
ww(n)=y(n);
ww(n-1)=y(n-1)-u(n-1,n)*ww(n);
for i=n-1:-1:2
    ww(i-1)=y(i-1)-u(i-1,i)*ww(i);
    ww(1)=y(1)-u(1,2)*ww(2);
end

% ww è la soluzione numerica 
w=linsolve(a,b);% risolvo anche il sistema con linsolv
plot(x,w,"xr",x, ww,"ob",x,z)
% vado a definire i coefficicenti p,q ed r 
%err=abs(z-ww)
err= abs((ww(5)-z(5)))
function z=fp(x) 
    z=-4*x^(-1);
end
function z=fq(x)
    z=-2*x^(-2);
end
function z=fr(x)
    z=2*x^(-2)*log(x);
end
% soluzione esatta 
function z=zsol(x)
    z=-2*x^(-2)+4*x^(-1)+log(x)-3/2;
end

