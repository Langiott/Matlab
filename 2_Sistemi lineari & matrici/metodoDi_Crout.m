clear all;
close all;
clc;
% costruisco la matrice
a = [10, -1, 0, 0; -2, 4, -2, 0; 0, -1, 2, -3; 0, 0, 2, -3]
% dichiaro i termini noti
d = [0; -1; 1.5; 2]
fprintf ("METODO CALCOLATO CON linsolve: \n");
linsolve(a,d)% ottengo la soluzione esatta del sistema
n=length(a);% la matrice è quadrata e ottengo la sua dimensione
l(1,1)=a(1,1);% costruisco la matrice diagonale
u(1,2)=a(1,2)/l(1,1);%  costruisco la matrice sovradiagonale
y(1)=d(1)/l(1,1);% faccio il rapporto tra il termine dono e i termini diagonali
for i=2:n
    l(i,i-1)=a(i,i-1);
end
for i=2:n-1
    l(i,i-1)=a(i,i-1);
    l(i,i)=a(i,i)-l(i,i-1)*u(i-1,i);
    u(i,i+1)=a(i,i+1)/l(i,i);
    y(i)=(d(i)-l(i,i-1)*y(i-1))/l(i,i);
end
l(n,n)=a(n,n)-l(n,n-1)*u(n-1,n);
y(n)=(d(n)-l(n,n-1)*y(n-1))/l(n,n);
x(n)=y(n);
x(n-1)=y(n-1)-u(n-1,n)*x(n);
for i=n-1:-1:2
    x(i-1)=y(i-1)-u(i-1,i)*x(i);
    x(1)=y(1)-u(1,2)*x(2);
end
fprintf ("METODO CALCOLATO CON METODO DI CROUT: \n");
for i=1:n
    fprintf ("%d\t%f\n", i, x(i));
end