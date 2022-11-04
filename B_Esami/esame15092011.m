clear all;
close all;
clc;
%% esercizio 1)

% determiniamo il numero di radici

sa=-pi/2;
sb= pi/2;
tol=1.e-2;
nr=100;
% grafico

x0=linspace(sa,sb,nr);
y0=f(x0);
%figure;
%plot(x0,y0);
fprintf("il sistema ammette una sola radice\n");

% definiamo lo zero con Nwton Raphson

% determiniamo la condizione iniziale 
xi=0.6;
zeros=NR(xi,tol,nr);
fzeros=f(zeros);
fprintf("lo zero si trova in %f\t\n " ,zeros);
fprintf("la funzione in  %f\t vale %d\t\n " ,zeros,fzeros);


%% esercizio 2)

% riportiamo i nodi inanzitutto
a=0;
b=3;
x= [0 0.1 0.2 0.3 2 2.5];
y= fg(x);
y2= [0 0.5 1 1.5 2 2.5];




xx=linspace(a,b,10000);% costruisco i vettore delle x della funzione interpolante
for i=1:10000
    %y1(i)=f(xx(i));
    yy(i)=interpol(x,y,xx(i));% salvo la funzione interpolata 
end



hold on;
plot(xx,yy);
plot(x,y,x,y2);





%% funzioni

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

function w=interpol(x,y,z)
%x,y vettori dei nodi
%z punto in cui interpolare
    sum=0;
    n=length(x);
    for k=1:n
        num=1;
        den=1;
        for j=1:k-1
            num=num*(z-x(j));
            den=den*(x(k)-x(j));
        end
        for j=k+1:n
            num=num*(z-x(j));
            den=den*(x(k)-x(j));
        end
        lnk=num/den;
        sum=sum+y(k)*lnk;
    end
    w=sum;
end

function y0=f(x)
    y0=exp(-x.^2)-sin(x);
end
% derivata 
function yp=fp(x)
    yp=-2*exp(-x.^2).*x - cos(x);
end

function yo=fg(x)
    yo=exp(-x.^2);
end
    