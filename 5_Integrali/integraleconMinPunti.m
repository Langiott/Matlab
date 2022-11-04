clear all;
close all;
clc;

%% estremi

a=0.0;
b=2*pi;
x=linspace(a,b,1001);
y=zeros(1,1001);

% costruisco la funzione 
for i=1:1001
    y(i)=f(x(i));
end
plot(x,y)
line([a,b],[0,0])
hold on

% costruisco queto vettore
p=0:pi/10:2*pi; % scelgo questo passo di ampionamento per valutare l'integrale
tol=1.0e-6;
err=1.0;
intold=0;
n=4; % punti iniziali
while err>tol
    ints=0;
    % costruisco i pesi
    ws(1)=1.0;
    ws(n+1)=1.0;
    for i=2:2:n-2
        ws(i)=4.0;
        ws(i+1)=2.0;
    end
    ws(n)=4.0;
    
   for k=1:20
       % valuto il valore dell' integrale nel seguente intervallino
        x1=p(k);
        x2=p(k+1);
        
        h=(x2-x1)/n;% passo di discretizzazione
        
        x=linspace(x1,x2,n+1);
        y=zeros(1,n+1);
        
        for i=1:n+1
            y(i)=f(x(i));
        end
        
        ints=ints+h*sum(ws.*y)/3.0
    end
    err=abs(ints-intold);
    fprintf("%d\t%f\t%e\n", n, ints, err);
    intold=ints;
    
    % che ho valutato il valore dell'integrale se soddisfa le condizioni
    % sull'errore
    n=n+2;
end

function z=f(x)
    z=exp(-x)*sin(10*x);
end


