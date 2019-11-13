%% NOMBRE: FLAMEO3DDOMINIOFRECUENCIA.m 
% CALCULO DE LA ECUACION AEROELASTICA Y CALCULO DE LAS CONDICIONES DE FLAMEO Y 
% VELOCIDAD DE DIVERGENCIA ESTATICA EN EL DOMINIO DE LA FRECUENCIA PARA UN ALA SEMIEMPOTRADA 

%% Parámetros de entrada 
% Cuerda del ala en voladizo 
c=0.2; 
% Envergadura del ala en voladizo 
s=0.6; 
% Espesor del ala en voladizo 
h=0.001; 
% Densidad del aluminio 
rhom=2770; 
% Módulo de Young del aluminio 
Eyoungs=68.7e9; 
% Coeficiente de Poisson del aluminio 
nu=0.33; 
% Rigidez del ala 
Dstiff=Eyoungs*h^3/12/(1-nu^2); 
% Densidad del aire 
rho=1.225; 
%% Malla de integración numérica 
% Número de puntos en la direccion de la cuerda y de la envergadura 
npoints=21; 
% Mallado en la direccion de la cuerda 
x=linspace(0,1,npoints); 
% Mallado en la direccion de la envergadura 
y=linspace(0,1,npoints); 
% Incremento en la direccion de la cuerda 
dx=x(2)-x(1); 
% Incremento en la direccion de la envergadura 
dy=y(2)-y(1); 
% Creacion de las matrices de malla 
[X,Y]=meshgrid(x,y); 
% La columna de X es igual a x y todas las columnas son identicas 
X=X'; 
% La fila de Y es igual a y y todas son identicas 
Y=Y'; 
% Seleccion de modos fuera del plano 
% Numero de modos en la direccion de la cuerda 
xmodos=4; 
% Numero de modos en la direccion de la envergadura 
ymodos=4; 
% Número de modos totales 
nmodos_out=xmodos*ymodos; 
% Calculo de las formas modales. 
[w,wxx,wyy,wxy,wx,wy]=formasmodales(xmodos,ymodos,nmodos_out,X,Y);
% No se han calculado los modos en el plano puesto que no aparecen 
% al linealizar las ecuaciones. 
% Dimensionamos las formas modales de las derivadas 
wx=wx/c; 
wy=wy/s; 
wxx=wxx/c^2; 
wyy=wyy/s^2; 
wxy=wxy/c/s; 
% Dimensionamos la malla numerica 
X=X*c; 
Y=Y*s; 
dx=dx*c; 
dy=dy*s; 
% Solo usaremos los modos fuera del plano 
nmodos=nmodos_out; 
% Establecemos las ecuaciones del movimiento 
% Establecer la matriz de masa estructural 
A=zeros(nmodos); 
for i=1:nmodos 
    for j=1:nmodos 
        % Calculo del integrando 
        integrand=rhom*h*squeeze(w(i,:,:).*w(j,:,:)); 
        % Integracion numerica (regla de Simpson) 
        A(i,j)=simpson2D(integrand,dx,dy); 
    end
end
% Establecer la matriz de rigidez estructural 
E=zeros(nmodos); 
for i=1:nmodos 
    for j=1:nmodos 
        % Calculo del integrando 
        integrand=Dstiff*squeeze(wxx(i,:,:).*wxx(j,:,:)+ wyy(i,:,:).*wyy(j,:,:)+nu*(wxx(i,:,:).*wyy(j,:,:)+ wyy(i,:,:).*wxx(j,:,:))+2*(1-nu)*wxy(i,:,:).*wxy(j,:,:)); 
        % Integracion numerica (regla de Simpson nuevamente) 
        E(i,j)=simpson2D(integrand,dx,dy); 
    end
end
%% Método del UVLM 
% Numero de paneles en la direccion de la cuerda 
mv=30; 
% Numero de paneles en la direccion de la envergadura 
nv=30; 
% Modelado de la aerodinamica del ala 
% El doble de la envergadura del ala en voladizo 
b=2*s; 
% Mitad de la cuerda del ala en voladizo 
cm=c/2; 
% Relacion de aspecto del ala aerodinamica 
AR=b/c; 
% Ángulo de ataque (Debe ser cero. No esta implementado el metodo para 
% otros angulos de ataque)
alpha=0*pi/180; 
% Longitud de la estela 
cw=10*c; 
% Multiplo de la longitud de la cuerda del ala 
% Calculo de los puntos de las esquinas de los paneles del ala 
% Cálculo de los nodos del ala 
Nodosala = MatrizNodos(mv,nv,0,0,b,c); 
% Coordendas de los nodos en la dirección de la envergadura 
nodosy=Nodosala(:,5)'; 
% Matriz de coordenadas de los nodos en la dirección de la envergadura 
yp=vec2mat(nodosy,nv+1); 
% Coordenadas de los nodos en la dirección de la cuerda 
nodosx=Nodosala(:,4)'; 
% Matriz de coordenadas de los nodos en la dirección de la cuerda 
xp = vec2mat(fliplr(nodosx),nv+1); 
% Coordenadas de los nodos en la dirección perpendicular al ala 
zp=zeros(mv+1,nv+1); 
% Calculo de las esquinas de anillos de torbellinos del ala 
xv = zeros(mv+1,nv+1); 
yv = zeros(mv+1,nv+1); 
zv = zeros(mv+1,nv+1); 
for j=1:nv+1 
    for i=1:mv 
        % Los anillos de torbellinos estan desplazados corriente abajo a un 
        % cuarto de su cuerda 
        xv(i,j)= xp(i,j)+(xp(i+1,j)-xp(i,j))/4; 
        yv(i,j)= yp(i,j); 
        zv(i,j)= zp(i,j)+(zp(i+1,j)-zp(i,j))/4; 
    end
    % Torbellino de salida de los anillos de torbellino del borde de fuga 
    xv(mv+1,j) = xp(mv+1,j)+ 0.25*c/mv; 
    yv(mv+1,j) = yp(mv+1,j); 
    zv(mv+1,j)= zp(mv+1,j); 
end
% Calculo de las esquinas de los anillos de torbellinos de la estela 
% Numero de anillos de torbellino en la direccion de la cuerda 
mw=cw*mv/c; 
% Calculo de los nodos de la estela 
Nodosestela = MatrizNodos(mw,nv,0,0,b,cw); 
% Posicionamiento de la estela en la geometria global 
wakepoints=min(xv(mv+1,:))+Nodosestela(:,4)'; 
% Coordenada x más cercana al ala 
xminestela=min(wakepoints); 
% Diferencia en x ala-estela 
difxalaestela=xminestela-min(xv(mv+1,:)); 
% Coordenadas en x de la estela unida al ala 
wakepoints=wakepoints-difxalaestela; 
% Coordenadas de los nodos, direccion de la cuerda 
xw = vec2mat(fliplr(wakepoints),nv+1); 
% Coordenadas de los nodos, direccion de la envergadura 
yw=repmat(yp(1,:),mw+1,1); 
% Coordenadas de los nodos, direccion perpendicular al plano de la estela 
zw=zeros(mw+1,nv+1); 
%% Representación gráfica del ala aerodinamica y de la estela 
figure(100) 
hmesh1=mesh(xp,yp,zp,'Edgecolor','b','FaceColor','none');
hold on 
hmesh2=mesh(xv,yv,zv,'Edgecolor','g','FaceColor','none'); 
hmesh3=mesh(xw(1:2*mv,:),yw(1:2*mv,:),zw(1:2*mv,:),'Edgecolor','r','FaceColor','none'); 
xlabel('x') 
ylabel('y') 
zlabel('z') 
ylim([-b/2 b/2]) 
legend([hmesh1(1) hmesh2(1) hmesh3(1)],'Paneles Geometricos','Anillos del ala','Anillos de la estela') 
%% Calculo de los vectores normales 
% El ala es siempre plana. Vectores normales apuntan en la direccion z 
nx = zeros(mv,nv); 
ny = zeros(mv,nv); 
nz = ones(mv,nv); 
% Dimensionamos nz en una matriz (mv*nv,3) 
Nall=[reshape(nx',mv*nv,1) reshape(ny',mv*nv,1) reshape(nz',mv*nv,1)]; 
% Puntos de control y areas de los anillos de torbelinos 
% Inicicialización de vectores 
% Longitud de los anillos de torbellinos en la direccion de la cuerda 
Dx = zeros(mv,nv); 
% Longitud de los anillos de torbellinos en la direccion de la envergadura 
Dy = zeros(mv,nv); 
% Coordenadas de los puntos de control 
xc = zeros(mv,nv); 
yc = zeros(mv,nv); 
zc = zeros(mv,nv); 
for i=1:mv 
    for j=1:nv 
        % Calculamos las posiciones de los puntos de control 
        dxv = (xp(i+1,j) - xp(i,j)+xp(i+1,j+1) - xp(i,j+1))/2; 
        dzv = (zp(i+1,j) - zp(i,j)+zp(i+1,j+1) - zp(i,j+1))/2; 
        xc(i,j) = (xp(i,j)+xp(i+1,j)+xp(i,j+1)+xp(i+1,j+1))/4+dxv/4; 
        yc(i,j) = (yp(i,j)+yp(i+1,j)+yp(i,j+1)+yp(i+1,j+1))/4; 
        zc(i,j) = (zp(i,j)+zp(i+1,j)+zp(i,j+1)+zp(i+1,j+1))/4+dzv/4; 
        % Calculamos las longitudes de los torbellinos de frontera 
        % Dx y Dy deben ser simetricos respecto a la linea central del ala 
        if j <= nv/2 
            Dx(i,j)=sqrt((xv(i+1,j)-xv(i,j))^2+(yv(i+1,j)-yv(i,j))^2+(zv(i+1,j)-zv(i,j))^2); 
            Dy(i,j)=sqrt((xv(i,j+1)-xv(i,j))^2+(yv(i,j+1)-yv(i,j))^2+(zv(i,j+1)-zv(i,j))^2); 
        else
            Dx(i,j)=sqrt((xv(i+1,j+1)-xv(i,j+1))^2+(yv(i+1,j+1)-yv(i,j+1))^2+(zv(i+1,j+1)-zv(i,j+1))^2); 
            Dy(i,j)=sqrt((xv(i,j)-xv(i,j+1))^2+(yv(i,j+1)-yv(i,j))^2+(zv(i,j+1)-zv(i,j))^2); 
        end
    end
end
% Areas de los anillos de torbellinos 
Sarea=Dx.*Dy; 
% Redimensionamos todas las matrices mv x nv en vectores mv*nv x 1 
nvec=[reshape(nx',mv*nv,1) reshape(ny',mv*nv,1) reshape(nz',mv*nv,1)]; 
Dyvec=reshape(Dy',mv*nv,1); 
Svec=reshape(Sarea',mv*nv,1); 
% Calculamos la influencia de los anillos del ala en el ala 
disp('Calculando la matriz de coeficientes de influencia del ala') 
% Inicializamos la matriz de coeficientes de influencia del ala 
Ab=zeros(mv*nv,mv*nv); 
for ic=1:mv 
    for jc=1:nv 
        for i=1:mv 
            for j=1:nv 
                % Velocidad inducida por el primer segmento de torbellino 
                uvw1=segmentotorbellino(xc(ic,jc),yc(ic,jc),zc(ic,jc),xv(i,j),yv(i,j),zv(i,j),xv(i,j+1),yv(i,j+1),zv(i,j+1)); 
                % Velocidad inducida por el segundo segmento de torbellino 
                uvw2=segmentotorbellino(xc(ic,jc),yc(ic,jc),zc(ic,jc),xv(i,j+1),yv(i,j+1),zv(i,j+1),xv(i+1,j+1),yv(i+1,j+1),zv(i+1,j+1)); 
                % Velocidad inducida por el tercer segmento de torbellino 
                uvw3=segmentotorbellino(xc(ic,jc),yc(ic,jc),zc(ic,jc),xv(i+1,j+1),yv(i+1,j+1),zv(i+1,j+1),xv(i+1,j),yv(i+1,j),zv(i+1,j)); 
                % Velocidad inducida por el cuarto segmento de torbellino 
                uvw4=segmentotorbellino(xc(ic,jc),yc(ic,jc),zc(ic,jc),xv(i+1,j),yv(i+1,j),zv(i+1,j),xv(i,j),yv(i,j),zv(i,j)); 
                % Suma de todas las velocidades inducidas 
                uvw=uvw1+uvw2+uvw3+uvw4; 
                % Calculamos los coeficientes de influencia 
                Ab((ic-1)*nv+jc,(i-1)*nv+j)=dot(uvw,nvec((ic-1)*nv+jc,:)); 
            end
        end
    end
end
%% Calculo de la influencia de los anillos de la estela en el ala 
disp('Calculando la matriz de coeficientes de influencia de la estela') 
% Inicializamos la matriz de coeficientes de influencia de la estela 
% Influencia en la direccion normal 
Aw = zeros(mv*nv,mw*nv); 
for ic=1:mv 
    for jc=1:nv 
        for i=1:mw 
            for j=1:nv 
                % Velocidad inducida por el primer segmento de torbellino 
                uvw1=segmentotorbellino(xc(ic,jc),yc(ic,jc),zc(ic,jc),xw(i,j),yw(i,j),zw(i,j),xw(i,j+1),yw(i,j+1),zw(i,j+1)); 
                % Velocidad inducida por el segundo segmento de torbellino 
                uvw2=segmentotorbellino(xc(ic,jc),yc(ic,jc),zc(ic,jc),xw(i,j+1),yw(i,j+1),zw(i,j+1),xw(i+1,j+1),yw(i+1,j+1),zw(i+1,j+1)); 
                % Velocidad inducida por el tercer segmento de torbellino 
                uvw3=segmentotorbellino(xc(ic,jc),yc(ic,jc),zc(ic,jc),xw(i+1,j+1),yw(i+1,j+1),zw(i+1,j+1),xw(i+1,j),yw(i+1,j),zw(i+1,j));
                % Velocidad inducida por el cuarto segmento de torbellino 
                uvw4=segmentotorbellino(xc(ic,jc),yc(ic,jc),zc(ic,jc),xw(i+1,j),yw(i+1,j),zw(i+1,j),xw(i,j),yw(i,j),zw(i,j)); 
                % Suma de todas las velocidades inducidas 
                uvw=uvw1+uvw2+uvw3+uvw4; 
                % Calculo del coeficiente de influencia 
                Aw((ic-1)*nv+jc,(i-1)*nv+j)=dot(uvw,nvec((ic-1)*nv+jc,:)); 
            end
        end
    end
end
%% Matrices de propagacion de la estela 
Pb=[zeros(nv,(mv-1)*nv) eye(nv,nv);zeros((mw-1)*nv,(mv-1)*nv) zeros((mw-1)*nv,nv)]; 
Pw=[zeros(nv,(mw-1)*nv) zeros(nv,nv);eye((mw-1)*nv,(mw-1)*nv) zeros((mw-1)*nv,nv)]; 
Pb=sparse(Pb); 
Pw=sparse(Pw); 
% Matrices de coeficientes para el calculo de la sustentacion 
Gy=eye(mv*nv)+[zeros(nv,(mv-1)*nv) zeros(nv,nv);-eye((mv-1)*nv, (mv-1)*nv) zeros((mv-1)*nv,nv)]; 
Gy=sparse(Gy.*repmat(Dyvec,1,mv*nv)); 
GS=sparse(eye(mv*nv).*repmat(Svec,1,mv*nv)); 
% Calculo de las matrices modales en los puntos de control 
[wc,~,~,~,wcx]=formasmodales(xmodos,ymodos,nmodos,xc(:,nv/2+1:nv)/c, yc(:,nv/2+1:nv)/s); 
% Dimensionalizar las derivadas de las formas modales 
wcx=wcx/c; 
% Reflejar las formas modales (tecnica espejo) 
% Inicializar los vectores de las formas modales de ambos lados 
w2=zeros(nmodos,mv,nv); 
wx2=zeros(nmodos,mv,nv); 
for i=1:nmodos 
    w2(i,:,:)=[fliplr(squeeze(wc(i,:,:))) squeeze(wc(i,:,:))]; 
    wx2(i,:,:)=[fliplr(squeeze(wcx(i,:,:))) squeeze(wcx(i,:,:))]; 
end
%% Representacion de todas las formas modales determinadas 
for i=1:nmodos 
    figure(i) 
    hplot=mesh(xc,yc,squeeze(w2(i,:,:))); 
    set(hplot,'EdgeColor',[ 0 0 1]) 
    xlabel('x') 
    ylabel('y') 
    zlabel(['w_' num2str(i) '(x,y)']) 
    xlim([-0.6 0.6]) 
    ylim([-0.6 0.6]) 
    zlim([-20 20]) 
    view(-30,42) 
end
%% Matrices de transformacion modal 
% Redimensionar w2 y wx2 en matrices de estructura (mv*nv,nmodos) 
wall=zeros(mv*nv,nmodos);
wxall=zeros(mv*nv,nmodos); 
for i=1:nmodos 
    wall(:,i)=reshape(squeeze(w2(i,:,:))',mv*nv,1); 
    wxall(:,i)=reshape(squeeze(wx2(i,:,:))',mv*nv,1); 
end
%% Calculo de las matrices PI y Pc 
% Matric Pe(0) 
PI=repmat(eye(nv),mw,1); 
Pc=sparse([zeros(nv,(mv-1)*nv) eye(nv)]); 
% Cálculo de L0 
L0=(Gy*cos(alpha))*((Ab+Aw*PI*Pc)\Nall(:,3)); 
% Indices del elemento de L0 y pared perteneciente a la semiala derecha 
halfind=zeros(1,mv*nv/2); 
for i=1:mv 
    halfind((i-1)*nv/2+1:i*nv/2)=(i-1)*nv+(nv/2+1:nv); 
end
% Cálculo de Q0 
Q0=(L0(halfind).'*wall(halfind,:)).'; 
%% Cálculo de la velocidad de flameo 
% Incremento para el cálculo numérico del Jacobiano 
dU=1e-8; 
% Hipótesis inicial de la velocidad de flameo 
U=20; 
% Hipótesis inicial de la frecuencia reducidad de flameo 
k=0.2; 
%% Cálculo de la velocidad de flameo utilizando el método de Newton-Raphson 
cond=0; 
while cond == 0 
    % Cálculo del determinante del flameo 
    [~,FF]=MFAG(k,mv,nv,mw,cm,Gy,GS,Ab,Aw,Pc,wxall,wall,A,E,rho,U,halfind,alpha,Q0); 
    % Evaluación numérica del Jacobiano 
    Jac=zeros(2,2); 
    % Con respecto a U 
    U=U+dU; 
    [~,FF1]=MFAG(k,mv,nv,mw,cm,Gy,GS,Ab,Aw,Pc,wxall,wall,A,E,rho,U,halfind,alpha,Q0); 
    Jac(:,1)=(FF1-FF)/dU; 
    U=U-dU; 
    % Con respecto a k 
    k=k+dU; 
    [~,FF1]=MFAG(k,mv,nv,mw,cm,Gy,GS,Ab,Aw,Pc,wxall,wall,A,E,rho,U,halfind,alpha,Q0); 
    Jac(:,2)=(FF1-FF)/dU; 
    k=k-dU; 
    % Resolución del sistema 
    solma=-Jac\FF; 
    % Cálculo del criterio de convergencia 
    crit=sqrt(solma'*solma); 
    % Actualización de la velocidad aerodinámica 
    U=U+solma(1); 
    % Actualización de la frecuencia reducida 
    k=k+solma(2); 
    % Prueba de convergencia
    if crit < 1e-12 
        cond=1; 
    end
end
% Velocidad aerodinámica de flameo 
Uflut=U; 
% Frecuencia reducida de flameo 
kflut=k; 
% Frecuencia de flameo en rad/s 
omegaflut=k*U/cm; 
disp('Velocidad de Flameo (m/s)') 
disp(Uflut) 
disp('Frecuencia de Flameo(Hz)') 
disp(omegaflut/2/pi) 
% Hipótesis inicial para la velocidad aerodinámica de divergencia 
U=20; 
%% Cálculo de la matriz de fuerza aerodinámica generalizada, frecuencia cero 
Q10=MFAG(0,mv,nv,mw,cm,Gy,GS,Ab,Aw,Pc,wxall,wall,A,E,rho,U,halfind,alpha,Q0); 
% Calculo de Q0S para el determinante 
if alpha~=0 
    % Inicializamos la matriz Q0S 
    Q0S=zeros(length(Q0),length(Q0)); 
    % Rellenamos la matriz Q0S con las columnas de Q0 
    for i=1:length(Q0) 
        Q0S(:,i)=Q0; 
    end
end
%% Cálculo de la velocidad de divergencia utilizando el método de 
% Newton-Raphson 
cond=0; 
while cond == 0 
    if alpha==0 
        % Cálculo de la función objetivo 
        DD=det(E-rho*U^2*Q10(:,:,1)); 
        % Evaluación del Jacobiano numéricamente 
        U=U+dU; 
        DD1=det(E-rho*U^2*Q10(:,:,1)); 
        Jac=(DD1-DD)/dU; 
        U=U-dU; 
        % Resolución del sistema de Newton-Raphson 
        solma=-Jac\DD; 
        % Cálculo del criterio de convergencia 
        crit=sqrt(solma'*solma); 
        % Actualización de la velocidad aerodinámica 
        U=U+solma; 
        % Prueba de convergencia 
        if crit < 1e-12 
            cond=1; 
        end
    else
        % Cálculo de la función objetivo 
        DD=det((E-rho*U^2*Q10(:,:,1))./(-rho*U^2*Q0S*sin(alpha)))-1; 
        % Evaluación del Jacobiano numéricamente 
        U=U+dU; 
        DD1=det((E-rho*U^2*Q10(:,:,1))./(-rho*U^2*Q0S*sin(alpha)))-1;
        Jac=(DD1-DD)/dU; 
        U=U-dU; 
        % Resolución del sistema de Newton-Raphson 
        solma=-Jac\DD; 
        % Cálculo del criterio de convergencia 
        crit=sqrt(solma'*solma); 
        % Actualización de la velocidad aerodinámica 
        U=U+solma; 
        % Prueba de convergencia 
        if crit < 1e-12 
            cond=1; 
        end
    end
end
%% Velocidad aerodinámica de divergencia estática 
Udiv=U; 
disp('Velocidad de divergencia estatica (m/s)') 
disp(Udiv)

%%
function [Q1,FF]=MFAG(k,mv,nv,mw,cm,Gy,GS,Ab,Aw,Pc,wxall,wall,A,E,rho,U,halfind,~,~) 
% Cálculo de la matriz de fuerzas aerodinámicas generalizada para un ala rectangular usando el metodo del Vortex Lattice. 
% Tambien calcula el determinante del flameo. 
% Angulo de ataque debe ser cero alpha=0; 
% Calculo de la matriz Pe 
Pe=sparse(kron(exp(-1j*2*k/mv*(1:mw)'),eye(nv))); 
% Calculo de la matriz L1 
L1=(Gy*cos(alpha)+1j*k/cm*GS)*((Ab+Aw*Pe*Pc)\(wxall+1j*k/cm*wall)); 
% Calculo de la matriz Q1 
Q1=(L1(halfind,:).'*wall(halfind,:)).'; 
% Calculo de la matriz completa aeroelástica 
aeromat=-(k*U/cm)^2*A+E-rho*U^2*Q1; 
% Calculo del determinante de la matriz aeroelástica 
flut_det=det(aeromat); 
% Calculo de la función objetivo 
FF=[real(flut_det);imag(flut_det)]; 
end
