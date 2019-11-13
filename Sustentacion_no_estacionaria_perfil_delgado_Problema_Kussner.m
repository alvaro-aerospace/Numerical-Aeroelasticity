%% SOLUCIÓN NUMÉRICA NO ESTACIONARIA PARA PERFILES DELGADOS. 
%% PERFIL ENTRANDO EN UNA RAFAGA. FUNCION DE KUSSNER 
% 
% Autor: Álvaro Fernández Villar
% Este código forma parte del trabajo final de Máster: 
% Solución numérica de problemas aerodinámicos no estacionarios mediante el método de la malla de torbellinos.
% Consultar trabajo para más referencias.
%
%-----------------------------------------------------------CÓDIGO----------------------------------------------------------------------%

%% Datos de entrada 
% Numero de paneles en el perfil 
n=100; 
% Numero de paneles en la estela 
m=150; 
% Tiempo de calculo adimensional 
taumax=50; 
% Velocidad de rafaga vertical 
w0=1; 

%% Calculo de algunos valores necesarios 
% Longitud de los paneles en la estela 
dxw=taumax/m; 

%% Mallado del perfil y de la estela 
% Inicializacion de algunas matrices 
% Matriz de puntos de control 
xc = zeros(n,1); 
% Matriz de puntos de torbellino 
xg = zeros(n,1); 
% Vector de tiempos de simulacion 
tau=zeros(1,m); 
% Puntos extremos de los paneles 
xk=transpose(linspace(-1,1,n+1)); 
% Puntos de control de cada panel 
for j=1:n 
    xc(j) = 1/4 * xk(j)+3/4 *xk(j+1); 
end
% Puntos de torbellino de cada panel 
for j=1:n 
    xg(j) = 1/4 * xk(j+1)+3/4 *xk(j); 
end
% Puntos de torbellino en la estela 
j=0:m-1; 
xkw=ones(1,m)+dxw/4; 
xgw=xkw(1:end)+j*dxw; 
xgw=transpose(xgw); 
% Instantes de simulacion
dtau=dxw; 
j=0:m-1; 
%Vector de instantes de simulacion 
tau=tau+j*dtau; 

%% Matrices del sistema de ecuaciones 
% Inicializacion de matrices 
% Matriz de coeficientes aerodinamicos 
A = zeros(n,1); 
% Matriz aw de coeficientes de la estela 
aw=zeros(n,m); 
% Vector velocidad de la condicion de frontera 
wa=zeros(n,m); 
% Matriz aw (Coeficientes aerodinamicos de la estela para un instante j) 
for i=1:n 
    for j=1:m 
        aw(i,j)=(-1/(2*pi))*(1/(xc(i)-xgw(j))); 
    end
end
% Matriz A (Coeficientes aerodinamicos del perfil) 
for i=1:n 
    for j=1:n 
        A(i,j)=(-1/(2*pi))*(1/(xc(i)-xg(j))); 
    end
end
% Calculo del vector wa 
for t=1:n 
    for q=1:m 
        wa(t,q)=-w0*heaviside((tau(t)-xc(t)-1)); 
    end
end

%% Resolucion del sistema y calculo de los torbellinos Gamma
% Inicializacion de matrices 
r=ones(n,1); 
% Vector de torbellinos de la estela 
Gammaw=zeros(1,m); 
% Matriz de gradiente numerico de torbellinos 
dGamma=zeros(n,m); 
% Matriz intermedia para el calculo de las presiones 
T=zeros(n,n); 
% Matriz de coeficiente de presiones 
dcp=zeros(n,m); 
% Vector intermedio calculo del cl 
cluno=zeros(1,m); 
% Vector intermedio calculo del cl 
cldos=zeros(1,m); 
% Vector del coeficiente de sustentacion teorico 
clexacto=zeros(1,m); 
% Vector solucion del coeficiente de sustentacion cuasi-estacionario 
clquasi=zeros(1,m); 
% Primer instante de tiempo 
Gamma(:,1)=(A-aw(:,1)*r')\wa(:,1); 
Gammaw(1)=-r'*Gamma(:,1); % Resto de instantes de tiempo
for j=2:m 
    suma=0; 
    for k=1:j-1 
        suma=suma+(aw(:,1)-aw(:,j-k+1))*Gammaw(k); 
    end
    Gamma(:,j)=(A-aw(:,1)*r')\(wa(:,j)+suma); 
    Gammaw(j)=-r'*Gamma(:,j); 
    for k=1:j-1 
        Gammaw(j)=Gammaw(j)-Gammaw(k); 
    end
end
% Calculo del gradiente numerico de Gamma 
for k=1:n 
    dGamma(k,:)=gradient(Gamma(k,:))/dtau; 
end

%% Calculo del coeficiente de presiones y coeficientes aerodinamicos
% Matriz T para el calculo del coeficiente de presiones 
for i=1:n 
    for j=1:n 
        if i==j 
            T(i,j)=3/2; 
        elseif i<j 
            T(i,j)=0; 
        elseif i>j 
            T(i,j)=1; 
        end
    end
end
% Coeficiente de presiones dcp 
for j=1:m 
    dcp(:,j)=n*Gamma(:,j)+T*dGamma(:,j); 
end
% Coeficiente de sustentacion numerico cl 
for i=1:m 
    % Primera parte 
    cluno(i)=sum(Gamma(:,i)); 
    % Segunda parte 
    cldos(i)=sum(dGamma(:,i).*(1-xg)); 
end
% Coeficiente de sustentacion numerico 
cl=cluno+cldos; 
% Coeficiente de sustentacion dado por Kussner y cuasi-estacionario 
Chi=funcionkussner(tau); 
%% Calculo de la funcion de Kussner 
% Coeficiente analitico 
for j=1:m 
    clexacto(j)=(2*pi()*w0*Chi(j)); 
end
%% Representacion de resultados 
figure 
hold on; 
% Coeficiente de sustentacion obtenido por el metodo de los paneles 
plot(tau,real(cl),'Linewidth',2,'color','red')
% Coeficiente de sustentacion obtenido por la funcion de Kussner 
plot(tau,real(clexacto),'-d','MarkerIndices',1:10:length(clexacto), 'Linewidth',2,'color','black') 
title(['Coeficiente de Sustentación CL(\tau), n=' num2str(n) ', m=' num2str(m) ', w0=' num2str(w0) ', \Delta\xi=' num2str(dxw)]) 
xlabel('Tiempo adimensional \tau, tau=t/Tr, Tr=b/U\infty');
ylabel('Coeficiente de Sustentación CL(\tau)'); 
legend('Solucion numerica', 'Solucion Kussner');

%% Funcion de Kusnner
function [Chi] = funcionkussner(tau) 
%Calculo de la funcion de Kussner 
% Inicializacion de matrices 
Chi=zeros(1,length(tau)); 
% Funcion de Kussner 
for m=1:length(tau) 
    Chi(m)=(tau(m).^2+tau(m))/(tau(m).^2+2.82*tau(m)+0.80); 
end
end
