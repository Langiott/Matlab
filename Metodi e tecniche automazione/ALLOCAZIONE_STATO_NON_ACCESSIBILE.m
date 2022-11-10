%PROBLEMA DI ASSEGNAZIONE DEGLI AUTOVALORI A CICLO CHIUSO
%nell'ipotesi stato non accessibile

n=3;  %numero di variabili di stato
m=2; %numero di ingressi
p=2; %numero di uscite

A=[1.2 0 0;
   -0.2 0.4 0;
    0 0 1];
B=[1 0;
   0 0;
   0 1];
C=[1 1 0;
   0 0 1];

D=zeros(p,m);

%verifica delle proprieta' strutturali
Mr=ctrb(A,B);
Mo=obsv(A,C);
nr=rank(Mr); %il sistema e' raggiungibile se e solo se  nr=n 
no=rank(Mo);  %il sistema e' osservabile se e solo se no=n 

%IL SISTEMA IN CATENA APERTA E' RAGGIUNGIBILE ED OSSERVABILE
%IL PROBLEMA DI ASSEGNAZIONE DEGLI AUTOVALORI A CICLO CHIUSO AMMETTE
%SOLUZIONE

%APPLICO IL PRINCIPIO DI SEPARAZIONE PER DETERMINARE I GUADAGNI K e L
K=place(A,B,[0.8 0.7 0.6]);
Lt=place(A',C',[0.2 0.3 0.4]);
L=Lt';

%rappresentazione in spazio di stato del sistema a ciclo chiuso
Ac=[A -B*K;L*C A-B*K-L*C];
eig(Ac)
Bc=[B;B];
Cc=[C zeros(p,n)];
Dc=zeros(p,m);

%1?simulazione 
N=100; %numero di campioni
x0=[1;1;1;0;0;0];
r=ones(N,m); %ingresso di riferimento
[Y,X]= dlsim(Ac,Bc,Cc,Dc,r,x0);

figure;
subplot(211),
plot(Y(:,1)),title('prima componente di y(k)');
grid
subplot(212)
plot(Y(:,2)),title('seconda componente di y(k)');
grid


%confronto  componente per componente stato vero e stato stimato
k=1:1:100;  %campioni
figure;
subplot(311)
plot(k,X(:,1),'r',k, X(:,n+1),'b');
grid
subplot(312)
plot(k,X(:,2),'r',k, X(:,n+2),'b');
grid
subplot(313)
plot(k,X(:,3),'r', k, X(:,n+3),'b');
grid