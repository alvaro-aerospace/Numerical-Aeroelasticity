%% SOLUCIÓN NUMÉRICA ESTACIONARIA PARA PERFILES DELGADOS. 
%% PERFIL SIN CURVATURA CON ANGULO DE ATAQUE 

%% Datos de entrada 
% Numero de paneles 
n=200; 
% Angulo de ataque en grados 
alpha0=10; 
% Velocidad del flujo 
U=1; 
% Cuerda del perfil 
c=1; 
% Densidad 
rho=1; 
% Presion dinamica 
Que=0.5*rho*U^2; 
%% Calculos intermedios 
% Angulo de ataque en radianes 
alpha=deg2rad(alpha0); 
% Semicuerda del perfil 
b=c/2; 
% Longitud del panel 
dx=c/n; 
%% Inicializacion de matrices 
xc=zeros(n,1); 
xg=zeros(n,1); 
A=zeros(n,n); 
dL=zeros(n,1); 
dcp=zeros(n,1); 
dM=zeros(n,1); 
cl=zeros(n,1); 
%% Calculo de los puntos caracteristicos del perfil 
% Puntos extremos de cada panel 
xj = transpose(linspace(-b,b,n+1)); 
% Puntos de control de cada panel 
for j=1:n 
    xc(j) = 0.25*xj(j)+0.75*xj(j+1); 
end
% Puntos de torbellino 
for j=1:n 
    xg(j) = 0.25*xj(j+1)+ 0.75*xj(j); 
end
%% Matriz A de coeficientes aerodinamicos 
for i=1:n 
    for j=1:n 
        A(i,j)=-1/(2*pi)*(1/(xc(i)-xg(j))); 
    end
end
%% Vector a 
a=-U*transpose(alpha*ones(1,n)); 
%% Resolucion de la circulacion g 
g=A\a; 
%% Calculo de las presiones y la sustentacion 
for i=1:n 
    % Calculo diferencia de sustentacion 
    dL(i)=rho*U*g(i); 
    % Calculo distribucion coeficiente sustentacion 
    cl(i)=dL(i)/(Que*c); 
    % Calculo coeficiente de presiones 
    dcp(i)=dL(i)/(dx*Que); 
    % Calculo momento aerodinamico con respecto al borde de ataque 
    dM(i)=-dL(i)*(xg(i)+b); 
end
%% Sustentacion total 
L=sum(dL); 
%% Coeficiente de sustentacion total 
cL=sum(cl); 
%% Coeficiente de momento aerodinamico con respecto al borde de ataque 
Cm0=sum(dM)/(Que*c^2); 
%% Ploteo del coeficiente de presiones 
plot(xg/c,dcp,'Linewidth',2,'color','blue') 
title({['Perfil delgado sin curvatura'],['Distribucion del coeficiente de Presiones cp(\eta), N=' num2str(n) ',\alpha=' num2str(alpha0) 'º']})
xlabel('\eta=X/c') 
ylabel('Coeficiente de Presiones \Deltacp(X)') 
legend('Solucion Numérica')
