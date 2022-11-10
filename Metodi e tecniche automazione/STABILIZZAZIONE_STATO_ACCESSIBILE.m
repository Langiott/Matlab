
%STABILITA' INTERNA  nell' ipotesi stato accessibile
%il problema ammette soluzione nel caso in cui (A,B) e' (almeno) stabilizzabile


n=3;  %numero di variabili di stato
m=1; %numero di ingressi
p=3; %numero di uscite
A=[2 0 -1;2 1 0;0 0 0.5];
eig(A)
B=[1; 0; 0];
C=eye(n);
D=zeros(p,m);

Mr=ctrb(A,B);
nr=rank(Mr)
%il sottospazio raggiungibile ha dimensione nr=2 
%il sottospazio non raggiungibile ha dimensione nnr=n-nr=1

nnr=n-nr;
%il Problema di stabilizzazione ammette soluzione nel caso in cui il sistema e'
%stabilizzabile ovvero il sottospazio non raggiungibile e' stabile.

%decomposizione rispetto alla raggiungibilita'
[a b c T k]=ctrbf(A,B,C)
%verifico se il sottospazio non raggiungibile e' stabile
Anr=a(1:nnr,1:nnr); %prende la sottomatrice estrapolando le righe dalla 1 alla nnr. Stessa cosa prendendo le colonne
eig(Anr) 
%il sottospazio non raggiungibile e' stabile. 

 
%progetto Kr che alloca gli autovalori alla coppia (Ar,Br) 
%(al sottospazio raggiungibile) 

Ar=a(nnr+1:n,nnr+1:n);
Br=b(nnr+1:n,:);
Kr=place(Ar,Br,[0.6 0.7]);


Khat=[zeros(m,nnr),Kr];
% dal  sistema trasformato occorre tornare al sistema originale
K=Khat*T;
eig(A-B*K)

%rappresentazione in spazio di stato del sistema a ciclo chiuso

Ac=A-B*K;
Bc=B;
Cc=C;
Dc=zeros(p,m);

% simulazione
N=100; %numero di campioni
x0=[2;0;1]; %condizione iniziale non nulla
r=zeros(N,m); %ingresso di riferimento nullo. Calcolo l'evoluzione libera del sistema a ciclo chiuso
[Y,X]= dlsim(Ac,Bc,Cc,Dc,r,x0);
figure;
subplot(311),
plot(Y(:,1)),title('prima componente di y(k)');
grid
subplot(312)
plot(Y(:,2)),title('seconda componente di y(k)');
grid
subplot(313)
plot(Y(:,3)),title('terza componente di y(k)');
grid


