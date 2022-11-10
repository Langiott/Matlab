%STABILITA INTERNA A CICLO CHIUSO 
%STATO NON ACCESSIBILE

n=3;  %numero di variabili di stato
m=2; %numero di ingressi
p=1; %numero di uscite

A=[1.2 0 0;0.2 1 0;0 0 0.8];
B=[1 0;0 0;0 1];
C=[0 1 0];


D=zeros(p,m);

%verifica delle proprieta' strutturali (stabilizzabilita' e rilevabilita')
Mr=ctrb(A,B);
Mo=obsv(A,C);
nr=rank(Mr); %il sistema e' raggiungibile  
no=rank(Mo); %il sistema non e' osservabile 
nno=n-no % dim sottospazio non osservabile
%decompongo rispetto all'osservabilita' per vedere se il sistema e'
%rilevabile
[a b c T k]=obsvf(A,B,C);
 ano=a(1,1);
 eig(ano) % il sistema e' rilevabile
 %estrapolare le matrici (Ao,Co) del sottospazio osservabile
 ao=a(2:3,2:3);
 co=c(:,2:3);
 lot=place(ao',co',[0.2 0.3]);
 lo=lot';
 Lhat=[zeros(nno,p);lo]
 L=inv(T)*Lhat;
 eig(A-L*C)
 
K=place(A,B,[0.6 0.5 0.4]);


Ac=[A -B*K;L*C A-B*K-L*C];
eig(Ac)
Bc=[B;B];
Cc=[C zeros(p,n)];
Dc=zeros(p,m);

%1 simulazione 
N=50; %numero di campioni
x0=[1;0.5;2;0;0;0];
r=zeros(N,m); %ingresso di riferimento
[Y,X]= dlsim(Ac,Bc,Cc,Dc,r,x0);

figure;
plot(Y)

