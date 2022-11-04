clear all;
close all;
clc;
tol=1.0e-6;
n=101;% imposto arbitrariamente un numero di iterazione
x=linspace(1,2*pi,n);% mi grafico la funzione e determino la condizione iniziale xi
plot(x,f(x));
xi=zeros(1,n);
xi(1)=3; % condizione iniziale

zr=NR(xi,tol,n)

function nraphs=NR(x,tol,n)
z(1)=0;
    for i=1:n
        dx=f(x(i))/fp(x(i));
        x(i+1)=x(i)-dx;
        z(i)= x(i+1);
        if(abs(dx)<=tol)
            nf=i;
            z(i)=x(i+1);
            break
        end
    
    end
    nraphs=z(end);
end
function y=f(x)
    y=sin(x);
end
function y=fp(x)% derivata di y
    y=cos(x);
end


    