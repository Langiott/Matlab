clear all;
close all;
clc;
% costruzione della funzione
a=0; b=1; %estremi dell' integrale
toll=1.e-12;
[err,iS]=sims(a,b,toll)

function [errore,Simpson]=sims(a,b,toll)
ap=0;
np=4;
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
    errore =er;
   
end

function  g=f(x)
    g=exp(-x);
end


