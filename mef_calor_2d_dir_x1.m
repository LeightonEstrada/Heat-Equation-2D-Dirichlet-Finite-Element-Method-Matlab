%ecuacion calor 2d elementos finitos 2d
%condiciones dirichlet
%tipo: http://es.wikipedia.org/wiki/Ecuaci%C3%B3n_del_calor
%leighton estrada rayme 04140403 unmsm fcm cc 14.5 
clear all;clc;
%-------------------------------------------------------------------------%
%AJUSTES INICIALES
%-------------------------------------------------------------------------%
tmax=15;%tiempo maximo
m=500;%numero particiones parte temporal
tt=tmax/m;%incremento tiempo
alfa=0.01;%dato inicial
%-------------------------------------------------------------------------%
%MALLADO TRIANGULAR
%-------------------------------------------------------------------------%
n=10;%numero particiones espaciales
h=1/n;ef=4*n^2;pu=(n+1)^2+n^2;%ef puntos %pu=nodos mallas
for j=1:n
    for i=(j-1)*n+1:j*n
        t(4*i-3,1)=i+j-1;t(4*i-2,1)=i+j-1;t(4*i-1,1)=i+j;t(4*i,1)=(n+1)^2+i;
        t(4*i-3,2)=(n+1)^2+i;t(4*i-2,2)=i+j;t(4*i-1,2)=n+1+i+j;t(4*i,2)=n+1+i+j;
        t(4*i-3,3)=n+i+j;t(4*i-2,3)=(n+1)^2+i;t(4*i-1,3)=(n+1)^2+i;t(4*i,3)=n+i+j;
    end
end
for j=1:n+1
    for i=(j-1)*(n+1)+1:j*(n+1)
        p(i,1)=(i-1-(j-1)*(n+1))*h;p(i,2)=(j-1)*h;
    end
end
for j=1:n
    for i=(j-1)*n+1:j*n
        pp(i,1)=(i-1-(j-1)*n)*h+(h/2);pp(i,2)=(j-1)*h+(h/2);
    end
end
p((n+1)^2+1:(n+1)^2+n^2,:)=pp;
%trimesh(t,p(:,1),p(:,2))
%title([num2str(ef),' elementos triangulares y ',num2str(pu),' puntos de malla'])
%-------------------------------------------------------------------------%
%ELEMENTOS FINITOS 2D
%-------------------------------------------------------------------------%
%CASO n=10
%ind=ones(pu,1);%ajustes frontera
%ind(1:11,1)=0;ind(111:121,1)=0;
%ind([12,23,34,45,56,67,78,89,100],1)=0;ind([22,33,44,55,66,77,88,99,110],1)=0;
%CASO n=3
%ind=ones(pu,1);%ajustes frontera
%ind(1:4,1)=0;ind(13:16,1)=0;ind([5,8,9,12],1)=0;
ind=((p(:,1)>0 & p(:,1)<1)&(p(:,2)>0 & p(:,2)<1));%obtiene puntos interiores caso general dom cuadrado
Np=size(p,1);N=sum(ind);% Np=# de nodos; N=# de nodos interiores
in=zeros(Np,1);in(ind)=(1:N)';
%for i=1:Np
%    if ind(i,1)==1 
%        in(i,1)=sum(ind([1:i],1));
%    end
%end  
K=zeros(N,N);%matriz de rigidez
M=zeros(N,N);%matriz de masa
for i=1:size(t,1) %no equivocarse aca
j=t(i,1);k=t(i,2);l=t(i,3);
vj=in(j);vk=in(k);vl=in(l);
J=[p(k,1)-p(j,1), p(l,1)-p(j,1); p(k,2)-p(j,2), p(l,2)-p(j,2)];
ar=abs(det(J))/2; ar1=abs(det(J))/24;CC=ar/12; Q=inv(J'*J);
if vj>0 
    K(vj,vj)=K(vj,vj)+ar*sum(sum(Q));
    M(vj,vj)=M(vj,vj)+ar1*2;
end
if vk>0
    K(vk,vk)=K(vk,vk)+ar*Q(1,1);
    M(vk,vk)=M(vk,vk)+ar1*2;
end
if vl>0
    K(vl,vl)=K(vl,vl)+ar*Q(2,2);
    M(vl,vl)=M(vl,vl)+ar1*2;
end
if vj*vk>0
    K(vj,vk)=K(vj,vk)-ar*sum(Q(:,1)); K(vk,vj)=K(vj,vk);
    M(vj,vk)=M(vj,vk)+ar1; M(vk,vj)=M(vj,vk);
end
if vj*vl>0
    K(vj,vl)=K(vj,vl)-ar*sum(Q(:,2)); K(vl,vj)=K(vj,vl);
    M(vj,vl)=M(vj,vl)+ar1; M(vl,vj)=M(vj,vl);
end
if vk*vl>0
    K(vk,vl)=K(vk,vl)+ar*Q(1,2); K(vl,vk)=K(vk,vl);
    M(vk,vl)=M(vk,vl)+ar1; M(vl,vk)=M(vk,vl);
end
end
%-------------------------------------------------------------------------%
%CONDICIONES INICIALES EXPERIMENTAL
%-------------------------------------------------------------------------%
U=zeros(pu,m+1);
U(26:30,1)=4;U(36:42,1)=4;U(48:53,1)=4;U(60:64,1)=4;U(70:75,1)=4;U(80:86,1)=4;U(92:96,1)=4;
U(121+23:121+28,1)=4;U(121+33:121+38,1)=4;U(121+44:121+48,1)=4;U(121+54:121+58,1)=4;
U(121+63:121+68,1)=4;U(121+73:121+78,1)=4;
U(121+14:121+17,1)=2;U(121+84:121+87,1)=2;
U([121+39,121+49,121+59,121+69],1)=2;
%-------------------------------------------------------------------------%
UU=U;Ne=Np-N;
indx=~((p(:,1)>0 & p(:,1)<1)&(p(:,2)>0 & p(:,2)<1));%obtiene puntos frontera caso general dom cuadrado
ix=(1:Np)';ii=ix(indx);
for i=1:Ne
    iix(i)=ii(Ne-i+1);
end
%ind(1:11,1)=0;ind(111:121,1)=0;
%ind([12,23,34,45,56,67,78,89,100],1)=0;ind([22,33,44,55,66,77,88,99,110],1)=0;
%for i=[[121:-1:110],100,99,89,88,78,77,67,66,56,55,45,44,34,33,23,22,12,[11:-1:1]]
for i=iix
    UU(i,:)=[];
end
%-------------------------------------------------------------------------%
%CRANK NICHOLSON
%-------------------------------------------------------------------------%
gamma=(-alfa/2)*tt;
for i=2:m+1
    UU(:,i)=(eye(N)-gamma*M^-1*K)^-1*(eye(N)+gamma*M^-1*K)*UU(:,i-1);       
end
%adicion a la matriz U los datos dirichlet
U(ind,:)=UU;%esto es importante
%-------------------------------------------------------------------------%
%ANIMACION 3D
%-------------------------------------------------------------------------%
ttt=0:tt:tmax;
for j=1:length(ttt)
    h=sprintf('%5.2f',ttt(j));
    set(gcf,'renderer','zbuffer');
set(gca,'nextplot','replacechildren');
caxis manual;caxis([min(min(U)) max(max(U))]);trisurf(t,p(:,1),p(:,2),U(:,j));title({['Ecuacion Calor 2D (Leighton Estrada R. CC FCM UNMSM)'];[num2str(ef),' elementos triangulares y ',num2str(pu),' nodos de malla'];['Tiempo t = ',num2str(h),' s']});xlabel('x');ylabel('y');zlabel('U(x,y;t)');axis([0 1 0 1 min(min(U)) max(max(U))]);
XYZ(j)=getframe;
end
%fin