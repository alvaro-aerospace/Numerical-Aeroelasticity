%% NOMBRE: FLAMEODOMINIODELTIEMPO_CURVASFLAMEO.m
% CALCULA LAS MATRICES LINEALES DE LAS ECUACIONES DEL MOVIMIENTO DE UN ALA
% SEMIEMPOTRADA EN FLUJO INCOMPRESIBLE Y LA FUNCION DE APROXIMACION RACIONAL
% RFA. DOMINIO DEL TIEMPO. CONDICIONES DE FLAMEO Y CURVAS DE FLAMEO

FLAMEO3DDOMINIOFRECUENCIA

%% Seleccion de los valores de frecuencias reducidas
kvec=[0 0.001 0.01 0.1 0.2 0.3 0.4 0.5 0.7 1.0]';
%% Creacion de la matriz de fuerza aerodinamica generalizada para todas las
% frecuencias reducidas
Qmat=zeros(nmodos,nmodos,length(kvec));
for ik=1:length(kvec)
k=kvec(ik);
Qmat(:,:,ik)=MFAG(k,mv,nv,mw,cm,Gy,GS,Ab,Aw,Pc,wxall,wall,A,E,rho,U,...
halfind,alpha,Q0);
end
%% Seleccion del numero de retardos aerodinamicos
nl=2;
%% Calculo de los coeficientes de los retardos
gammai=1.7*max(kvec)*((1:nl)/(nl+1)).^2;
%% Funcion de aproximacion racional de fuerza aerodinamica generalizada
[err,S,Qbar]=RFA(nmodos,nl,gammai,Qmat,kvec);
%% Numero de estados totales
nstates=(2+nl)*nmodos;
%% Calculo de la estabilidad del sistema lineal para varias velocidades
% aerodinamicas
% Seleccion del rango de velocidades aerodinamicas
Uvec=linspace(1,40,401);
% Numero de velocidades aerodinamicas
nU=length(Uvec);
% Inicializacion del vector de autovalores
eigvals=zeros((2+nl)*nmodos,nU);
% Inicializacion del vector de frecuencias reducidas
omegaRFA=zeros((2+nl)*nmodos,nU);
% Inicializacion del vector de relaciones de amortiguamiento
dampRFA=zeros((2+nl)*nmodos,nU);
for iU=1:nU
% Calculo de la matriz del sistema
H=FlameoRFA_linmat(Uvec(iU),rho,cm,A,E,S,Q0,gammai,nmodos,nl);
% Calculo de los autovalores del sistema
eigvals(:,iU)=eig(H);
% Ordenacion de los autovalores
[dum,I]=sort(imag(eigvals(:,iU)),'descend');
eigvals(:,iU)=eigvals(I,iU);
end
%%
figure(1)
hplot=plot(Uvec,abs(imag(eigvals)./2./pi),'k.',Uflut,omegaflut/2/pi,'mo',...
'linewidth',2,'markersize',4);
xlabel('U')
ylabel('\Im(\lambda)')
title('Diagrama U-\Im(\lambda)')
grid on;
ylim([0 100]) 
xlim(Uvec([1 end])) 
string3=['\omega_F= ',num2str(round(omegaflut/2/pi,2))]; 
legend(hplot([1 (nl+2)*nmodos+1]),{'\lambda',string3}) 

figure(2) 
hplot=plot(Uvec,real(eigvals),'k.',Uflut,0,'mo',Udiv,0,'bd','linewidth',2,'markersize',4); 
set(hplot((nl+2)*nmodos+2),'color',[0 0 1]) 
xlabel('U') 
ylabel('\Re(\lambda)') 
title('Diagrama U-\Re(\lambda)') 
grid on; 
ylim([-10 2]) 
xlim(Uvec([1 end])) 
string1=['U_F= ',num2str(round(Uflut,1))]; 
string2=['U_D= ',num2str(round(Udiv,1))]; 
legend(hplot([1 (nl+2)*nmodos+1 (nl+2)*nmodos+2]),{'\lambda',string1,string2}) 

figure(3) 
hplot=plot(Uvec,abs((eigvals)./2./pi),'k.',Uflut,omegaflut/2/pi,'mo','linewidth',2,'markersize',4); 
xlabel('U') 
ylabel('\omega_n(Hz)') 
title('Curva de flameo U- \omega_n') 
grid on; 
ylim([0 100]) 
xlim(Uvec([1 end])) 
string3=['\omega_F= ',num2str(round(omegaflut/2/pi,2))]; 
legend(hplot([1 (nl+2)*nmodos+1]),{'\lambda',string3}) 

figure(4) 
hplot=plot(Uvec,-real(eigvals)./(abs((eigvals)./2./pi)),'k.',Uflut,0,'mo',Udiv,0,'bd','linewidth',2,'markersize',4); 
set(hplot((nl+2)*nmodos+2),'color',[0 0 1]) 
xlabel('U') 
ylabel('\zeta') 
title('Curva de flameo U- \zeta') 
grid on; 
ylim([-2 2]) 
xlim(Uvec([1 end])) 
string1=['U_F= ',num2str(round(Uflut,1))]; 
string2=['U_D= ',num2str(round(Udiv,1))]; 
legend(hplot([1 (nl+2)*nmodos+1 (nl+2)*nmodos+2]),{'\lambda',string1,string2})

%%
function [err,S,Qbar]=RFA(nmodos,nl,gammai,Qmat,kvec) 
% Cálculo de la función racional de aproximacion de la matriz de fuerzas 
% aerodinámicas generalizadas. 
% Número de frecuencias reducidas 
nk=length(kvec);
% Inicialización de matrices de coeficientes para la funcion racional de 
% aproximación 
S=zeros(nmodos,nmodos,3+nl); 
% Coeficiente constante 
S(:,:,1)=Qmat(:,:,1); 
% Configuración del ajuste de la curva. El ajuste se lleva a cabo para cada 
% elemento de la matriz de fuerza aerodinámica generalizada por separado. 
% Creacion de la matriz del lado izquierdo. Esta matriz es la misma para 
% todos los elementos de la matriz de fuerza aerodinámica generalizada. 
LHS=zeros(nk,2+nl); 
LHS(:,1:2)=[1i*kvec -kvec.^2]; 
for inl=1:nl 
    LHS(:,2+inl)=1i*kvec./(1i*kvec+gammai(inl)); 
end
% Realizar el ajuste de la curva 
for i=1:nmodos 
    for j=1:nmodos 
        % Creacion el vector del lado derecho. Este vector tiene un valor 
        % diferente para cada elemento de la matriz de fuerza aerodinamica 
        % generalizada. 
        RHS=squeeze(Qmat(i,j,:))-S(i,j,1); 
        % Separacion de la parte real e imaginaria 
        solma=[real(LHS);imag(LHS)]\[real(RHS);imag(RHS)]; 
        % Asignar valores para los elementos de la matriz de coeficientes de la RFA 
        S(i,j,2)=solma(1); 
        S(i,j,3)=solma(2); 
        for inl=1:nl 
            S(i,j,3+inl)=solma(2+inl); 
        end
    end
end
% Cálculo de la matriz Q1 reconstruida 
Qbar=zeros(nmodos,nmodos,nk); 
for ik=1:nk 
    % Aproximación de Roger 
    Qbar(:,:,ik)=S(:,:,1)+S(:,:,2)*1i*kvec(ik)-S(:,:,3)*kvec(ik)^2; 
    for inl=1:nl 
        Qbar(:,:,ik)=Qbar(:,:,ik)+S(:,:,3+inl)*1i*kvec(ik)/(1i*kvec(ik)+gammai(inl)); 
    end
end
% Error de la curva de ajuste 
Qvec=reshape(Qmat,nmodos*nmodos*nk,1); 
Qbarvec=reshape(Qbar,nmodos*nmodos*nk,1); 
err=1/nk*sum(abs(Qvec-Qbarvec).^2./max([abs(Qvec).^2 ones(nmodos*nmodos*nk,1)],[],2));
end

%%
function [H,h0,Mbar]=FlameoRFA_linmat(U,rho,cm,A,E,S,Q0,gammai,nmodos,nl) 
% Calcula las matrices lineales de las ecuaciones del movimiento de primer orden 
% de un ala en voladizo en flujo incompresible con la función de aproximación 
% racional aerodinámica. 
% Matriz total de masa 
Mbar=A-rho*cm^2*S(:,:,3); 
% Matriz total de amortiguamiento 
Cbar=-rho*U*cm*S(:,:,2);
% Matriz total de rigidez 
Kbar=E-rho*U^2*S(:,:,1); 
% Matrices de estados aerodinámicos 
Sbar=zeros(nmodos,nl*nmodos); 
for inl=1:nl 
    Sbar(:,nmodos*(inl-1)+1:nmodos*inl)=-rho*U^2*S(:,:,3+inl); 
end
% Ecuación de estado aerodinámico 
WW0=zeros(nl*nmodos,2*nmodos); 
WW1=zeros(nl*nmodos,nl*nmodos); 
for inl=1:nl 
    WW0(nmodos*(inl-1)+1:nmodos*inl,1:nmodos)=eye(nmodos); 
    WW1(nmodos*(inl-1)+1:nmodos*inl,nmodos*(inl-1)+1:nmodos*inl)=-U*gammai(inl)/cm*eye(nmodos); 
end
% Matriz del sistema completo 
H=[-Mbar\Cbar -Mbar\Kbar -Mbar\Sbar; 
    eye(nmodos) zeros(nmodos,(nl+1)*nmodos); 
    WW0 WW1]; 
% Matriz de constantes 
h0=[-rho*U^2*(Mbar\Q0);zeros((1+nl)*nmodos,1)];
end
