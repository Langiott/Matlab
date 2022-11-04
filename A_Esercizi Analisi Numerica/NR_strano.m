close all;
clear all;
clc;
a=0;
b=2;

%%esercizio 1)

toll=1e-5;
x=linspace(0,2,100);
n=length(x);



% costruzione delle funzioni di cui vedere l'intersezione
y=f1(x);
y1=f2(x);
hold on
tt=y-y1;

toll=1e-5;
plot(x,y);
plot(x,y1);


%devo trovare l'intersezione quindi
x0=1.5;
xs(1)=x0;

z=0;
for i=1:n
    dx=f(xs(i))/fp(xs(i));
    xs(i+1)=xs(i)+dx;
    if(abs(dx)<=toll)
        nf=i;
        z(1)=xs(i+1);
        break
    end  
end

%% interfaccia 
fprintf ("il punto di intersezione tra le due funzioni vale %f:\n\n", z);

[err,sim]=sims(a,z,toll);
fprintf ("l'integrale calcolato tra %d e %f vale %f\ncon errore %d ",a,z,sim,err);



%% funzioni

function [errore,Simpson]=sims(a,b,toll)
ap=0;
np=14;
er=1;
    while er>toll
        h=(b-a)/np;
        x=linspace(a,b,np+1);
        for i=1:np+1
            y(i)=f(x(i));% creo i pu ti della funzione 
        end
        lung=length(y);
        % creo il vettore dei coefficienti per il metodo di simpson per i termini 
        % pari e dispari
        ws(1)=1.0;% inizializzo a 1 il primo termine dei coefficienti
        ws(np+1)=1.0;% inizializzo a 1 gli elementi successivi

        for i=2:2:np-2
            ws(i)=4.0;% inserisco 4 nei coefficienti parti
            ws(i+1)=2.0;% inserisco i 2 nei coefficienti pari
        end
        ws(np)=4;% questo è il penultimo elemento
        ints=h*sum(ws.*y)/3.0;% formula per simpson
        er=abs(ints-ap);
        ap=ints;
        np=np+2;
    end
    Simpson=ints;
    errore=er;
   
end


function h=f(x)
    h=f1(x)-f2(x);
end

% derivata della h funzione differenza
function h1=fp(x)
    h1=(4./(16*x.^2+1))-((sec(x)).^2)/4;
end

function y=f1(x)
    y=(1/4)*tan(x);
end

function y1=f2(x)
    y1=atan(4*x);
end