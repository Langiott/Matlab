clear all;
close all;
clc;

%% esercizio 1)
x12=0;
for j=1:16
    x12(j+1)=x12(j)+(1.25); % qui bisogna dichiarare gli indici successivi
end

% ho definito il vettore delle x 

% definisco il vettore delle y
y=[ 100 114  125 130 132 137 140  141  141 140  137  135 132 125 114 110 100];

% interpolo 

z=linspace(0,20,101);

for i=1:101
    intery(i)=interpol(x12,y,z(i));
end

hold on
plot (x12,y,"ro");
plot(z,intery);

%% esercizio 2)

a=[-4 1 -1 1;1 4 -1 -1;-1 -1 5 1;1 -1 1 3]
n=size(a,1)
b=[-2 -1 0 1]
ab=[a b'];
format long;
xtrue=linsolve(a,b')% soluzione quella vera
x=zeros(1,n);% inizializzo x con il vettore nullo zeros
err=1.0;
tol=1.0e-7;
k=1;
%applico il metodo ddi jacobi e ne verifico il raggio spettrale
while err>tol
%for k=1:nmax
    xold=x;% mi serve per confrontare i vecchi valori con qquelli anuovi per valutare l'errore
    for i=1:n
        sum=0;
        % mi vado a calcolare il membro di destra 
        % costruisco la sommatoria
        for j=1:i-1
            sum=sum+a(i,j)*xold(j);
        end
        % salto il termine i=j perchè è il membro di sinistra
        for j=i+1:n
            sum=sum+a(i,j)*xold(j);
        end
        % a questo punti aggiorno il valore diagonale 
        x(i)=(b(i)-sum)/a(i,i);
    end
    
    err=max(abs(x-xold));% valuto l'errore 
    %fprintf("%d\t%e\n", k, err);% mi stampo l'indice della iterazione 
    k=k+1;% incremento l'indice
end
% stampo le soluzioni ottenuta con il metodo vero,quella con il metodo
% linsolve e l'errore
for i=1:n
    fprintf("%f\t%f\t%e\n", x(i), xtrue(i), abs(x(i)-xtrue(i)));
end
fprintf ("\n Matrice T e il raggio spettrale \n\n");
d=zeros(n);
% controlliamo il raggio spettrale
for i=1:n
    d(i,i)=a(i,i);
end
r=a-d;
t=-inv(d)*r% costruisco la matrice t
fprintf (" \n\n");
lam=eig(t);% mi calcolo gli autovalori
rho=max(abs(lam))% e ne prendo il max
fprintf("il raggio spettrale è minore di 1 quindi\n");
fprintf("il metodo converge alla soluzione esatta\n");
% applico il metodo di gauss_saidel
% per sicurezza riassegno le variabili in modo da essere sicuri che il
% metodo funzioni
err=1.0;
tol=1.0e-7;% definisco la tolleranza 
x1=zeros(1,n); 
xoldn=x1;
k=1;

while err > tol
    xoldn=x1;
    for i=1:n
        sum=0;
        for j=1:i-1
            sum=sum+a(i,j)*x1(j);% sfrutto l'aggiornamento precedente
        end
        for j=i+1:n
            sum=sum+a(i,j)*xoldn(j);
        end
        x1(i)=(b(i)-sum)/a(i,i);
    end
    err=max(abs(x1-xoldn));
    %fprintf("%d\t%e\n", k, err);
    k=k+1;
    

end
for i=1:n
    fprintf("%f\t%f\t%e\n", x1(i), xtrue(i), abs(x1(i)-xtrue(i)));
end
% procedo con la costruzione della matrice T
fprintf ("\n Matrice t e raggio spettrale \n\n");
d=zeros(n);% inizializzo la matrice diagonale
for i=1:n
    d(i,i)=a(i,i);
end
l=zeros(n);% inizializzo la matrice triagolare inferiore
for i=1:n
    for j=1:i-1
        l(i,j)=-a(i,j);
    end
end
u=zeros(n);% Inizializzo la matrice triangolare superiore
for i=1:n
    for j=i+1:n
        u(i,j)=-a(i,j);
    end
end
t1=inv(d-l)*u
lam1=eig(t);
rho1=max(abs(lam1))
fprintf("il raggio spettrale è minore di 1 quindi il metodo converge");















%% funzioni
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


