% SOLUCION NUMERICA NO ESTACIONARIA PARA EL CASO TRIDIMENSIONAL. (UVLM) 
% CONDICIÓN DE FRONTERA: SALTO DE VELOCIDAD 

% Parámetros geométricos del ala 
% Dimensiones del ala 
% Longitud de la cuerda en el encastre 
cr=1; 
% Longitud de la envergadura alar 
b=3; 
% Flecha en el borde frontal (grados) 
LF=0; 
% Flecha en el borde posterior (grados) 
LR=0; 
% Salto de velocidad (condicion de frontera) 
Usalto=1; 
% Ángulo de ataque 
alpha=5*pi/180; 
% Parámetros de discretización de la geometría alar y de la estela 
% Número de paneles en la dirección de la envergadura 
nv=20; 
% Número de paneles en la dirección de la cuerda 
mv=16; 
% Densidad del aire 
rho=1.225; 
% Longitud de la estela en términos de cuerda en el encastre 
cw=10*cr; 
% Parámetros de discretización temporal 
% Número de pasos temporales 
ntimes=160; 
% Cálculo del paso temporal 
dt=cr/mv/(Usalto); 
% Creación del vector de tiempos 
ti=(0:ntimes-1)*dt; 
% Calculos de generación de la geometría global del problema 
% Cálculo de los nodos de los paneles del ala 
% Cálculo de los nodos del ala 
Nodosala = MatrizNodos(mv,nv,LF,LR,b,cr); 
% Coordendas de los nodos en la dirección de la envergadura 
nodosy=Nodosala(:,5)'; 
% Matriz de coordenadas de los nodos en la dirección de la envergadura 
yp=vec2mat(nodosy,nv+1); 
% Coordendas de los nodos en la dirección de la cuerda 
nodosx=Nodosala(:,4)'; 
% Matriz de coordenadas de los nodos en la dirección de la cuerda
xp = vec2mat(fliplr(nodosx),nv+1); 
% Coordenadas de los nodos en la dirección perpendicular al ala 
zp=zeros(mv+1,nv+1); 
% Cálculo de las esquinas de los anillos de torbellinos del ala 
xv = zeros(mv+1,nv+1); 
yv = zeros(mv+1,nv+1); 
zv = zeros(mv+1,nv+1); 
for j=1:nv+1 
    for i=1:mv 
        % Anillos de torbellino desplazados aguas abajo a 1/4 de su cuerda 
        xv(i,j)= xp(i,j)+(xp(i+1,j)-xp(i,j))/4; 
        yv(i,j)= yp(i,j); 
        zv(i,j)= zp(i,j)+(zp(i+1,j)-zp(i,j))/4; 
    end
    % Último segmento de los anillos de torbellino del borde de fuga 
    xv(mv+1,j) = xp(mv+1,j)+ 0.25*cr/mv; 
    yv(mv+1,j) = yp(mv+1,j); 
    zv(mv+1,j)= zp(mv+1,j); 
end
% Esquinas de los anillos de torbellino de la estela 
% Numero de anillos de torbellino en la direccion de la cuerda 
mw=cw*mv/cr; 
% Calculo de los nodos de la estela 
Nodosestela = MatrizNodos(mw,nv,LR,0,b,cw); 
% Posicionamiento de la estela en la geometria global 
wakepoints=min(xv(mv+1,:))+Nodosestela(:,4)'; 
% Coordenada x más cercana al ala 
xminestela=min(wakepoints); 
% Diferencia en x ala-estela 
difxalaestela=xminestela-min(xv(mv+1,:)); 
% Coordenadas en x de la estela unida al ala 
wakepoints=wakepoints-difxalaestela; 
% Matriz de coordenadas de los nodos en la direccion de la cuerda 
xw = vec2mat(fliplr(wakepoints),nv+1); 
% Matriz de coordenadas de los nodos en la direccion de la envergadura 
yw=repmat(yp(1,:),mw+1,1); 
% Matriz de coordenadas de los nodos, direccion normal al plano del ala 
zw=zeros(mw+1,nv+1); 
% Representación gráfica del ala y de la estela 
figure(100) 
hmesh1=mesh(xp,yp,zp,'Edgecolor','black','FaceColor','none'); 
hold on 
hmesh2=mesh(xv,yv,zv,'Edgecolor','blue','FaceColor','none'); 
hmesh3=mesh(xw(1:mw,:),yw(1:mw,:), zw(1:mw,:),'Edgecolor','green','FaceColor','none'); 
xlabel('x') 
ylabel('y') 
zlabel('z') 
ylim([-b/2 b/2]) 
legend([hmesh1(1) hmesh2(1) hmesh3(1)], 'Paneles Geometricos','Anillos del ala','Anillos de la estela') 
pbaspect([1 b/cr 1]) 
% Puntos de control, dimensiones y áreas de los anillos de torbellino 
% Inicialización de las longitudes de los anillos. Direccion cuerda 
Dx = zeros(mv,nv);
% Inicialización de las longitudes de los anillos. Direccion envergadura 
Dy = zeros(mv,nv); 
% Inicialización de los vectores de coordenadas de los puntos de control 
% en los anillos de torbellino del ala 
xc = zeros(mv,nv); 
yc = zeros(mv,nv); 
zc = zeros(mv,nv); 
% Inicializacion de los vectores de coordenadas de los puntos de control 
% en los anillos de torbellino de la estela 
xwc = zeros(mw,nv); 
ywc = zeros(mw,nv); 
zwc = zeros(mw,nv); 
for i=1:mv 
    for j=1:nv 
        % Cálculo de las posiciones de los puntos de control 
        dxv = (xp(i+1,j) - xp(i,j)+xp(i+1,j+1) - xp(i,j+1))/2; 
        dzv = (zp(i+1,j) - zp(i,j)+zp(i+1,j+1) - zp(i,j+1))/2; 
        xc(i,j) = (xp(i,j)+xp(i+1,j)+xp(i,j+1)+xp(i+1,j+1))/4+dxv/4; 
        yc(i,j) = (yp(i,j)+yp(i+1,j)+yp(i,j+1)+yp(i+1,j+1))/4; 
        zc(i,j) = (zp(i,j)+zp(i+1,j)+zp(i,j+1)+zp(i+1,j+1))/4+dzv/4; 
        % Cálculo de las longitudes de los torbellinos 
        % Asegurarse que Dx y Dy quedan simétricos alrededor del encastre 
        if j <= nv/2 
            Dx(i,j)=sqrt((xv(i+1,j)-xv(i,j))^2+(yv(i+1,j)-yv(i,j))^2+(zv(i+1,j)-zv(i,j))^2); 
            Dy(i,j)=sqrt((xv(i,j+1)-xv(i,j))^2+(yv(i,j+1)-yv(i,j))^2+(zv(i,j+1)-zv(i,j))^2); 
        else
            Dx(i,j)=sqrt((xv(i+1,j+1)-xv(i,j+1))^2+(yv(i+1,j+1)-yv(i,j+1))^2+(zv(i+1,j+1)-zv(i,j+1))^2); 
            Dy(i,j)=sqrt((xv(i,j)-xv(i,j+1))^2+(yv(i,j+1)-yv(i,j))^2+(zv(i,j+1)-zv(i,j))^2); 
        end
    end
end
% Calculo de los puntos de control en la estela 
for i=1:mw 
    for j=1:nv 
        % Cálculo de las posiciones de los puntos de control 
        xwc(i,j) = (xw(i,j)+xw(i+1,j))/2; 
        ywc(i,j) = (yw(i,j)+yw(i,j+1))/2; 
        zwc(i,j) = zw(i,j); 
    end
end
% Cálculo de los vectores normales 
nx = zeros(mv,nv); 
ny = zeros(mv,nv); 
nz = ones(mv,nv);
% Areas de los anillos de torbellinos 
S=Dx.*Dy; 
% Conversión de las matrices mv x nv a vectores mv*nv x 1 
nvec=[reshape(nx',mv*nv,1) reshape(ny',mv*nv,1) reshape(nz',mv*nv,1)]; 
Dyvec=reshape(Dy',mv*nv,1); Svec=reshape(S',mv*nv,1); 
% Calculo de la influencia de los torbellinos del ala en el ala 
disp('Calculando la matriz de coeficientes de influencia del ala') 
% Inicializando matriz de coeficientes de influencia del ala 
Ab=zeros(mv*nv,mv*nv); 
for ic=1:mv 
    for jc=1:nv 
        for i=1:mv 
            for j=1:nv 
                % Velocidad inducida por el primer segmento de torbellino 
                uvw1=segmentotorbellino(xc(ic,jc),yc(ic,jc),zc(ic,jc),xv(i,j),yv(i,j),zv(i,j),xv(i,j+1),yv(i,j+1),zv(i,j+1)); 
                % Velocidad inducida del segundo segmento de torbellino 
                uvw2=segmentotorbellino(xc(ic,jc),yc(ic,jc),zc(ic,jc),xv(i,j+1),yv(i,j+1),zv(i,j+1),xv(i+1,j+1),yv(i+1,j+1),zv(i+1,j+1)); 
                % Velocidad inducida por el tercer segmento de torbellino 
                uvw3=segmentotorbellino(xc(ic,jc),yc(ic,jc),zc(ic,jc),xv(i+1,j+1),yv(i+1,j+1),zv(i+1,j+1),xv(i+1,j),yv(i+1,j),zv(i+1,j)); 
                % Velocidad inducida por el cuarto segmento de torbellino 
                uvw4=segmentotorbellino(xc(ic,jc),yc(ic,jc),zc(ic,jc),xv(i+1,j),yv(i+1,j),zv(i+1,j),xv(i,j),yv(i,j),zv(i,j)); 
                % Suma de todas las velocidades inducidas 
                uvw=uvw1+uvw2+uvw3+uvw4; 
                % Cálculo de los coeficientes de influencia 
                Ab((ic-1)*nv+jc,(i-1)*nv+j)=dot(uvw,nvec((ic-1)*nv+jc,:)); 
            end
        end
    end
end
% Cálculo de la influencia de la estela en el ala 
disp('Calculando la matriz de coeficientes de influencia de la estela') 
% Inicialización de la matriz de coeficientes de influencia de la estela 
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
                % Cálculo de los coeficientes de influencia 
                Aw((ic-1)*nv+jc,(i-1)*nv+j)=dot(uvw,nvec((ic-1)*nv+jc,:)); 
            end
        end
    end
end
% Matrices de propagación de la estela 
Pb=[zeros(nv,(mv-1)*nv) eye(nv,nv);zeros((mw-1)*nv,(mv-1)*nv) zeros((mw-1)*nv,nv)]; 
Pw=[zeros(nv,(mw-1)*nv) zeros(nv,nv);eye((mw-1)*nv,(mw-1)*nv) zeros((mw-1)*nv,nv)]; 
Pb=sparse(Pb); 
Pw=sparse(Pw); 
% Matrices de coeficientes para el calculo de la sustentacion 
Gy=eye(mv*nv)+[zeros(nv,(mv-1)*nv) zeros(nv,nv);-eye((mv-1)*nv, (mv-1)*nv) zeros((mv-1)*nv,nv)]; 
Gy=sparse(Gy.*repmat(Dyvec,1,mv*nv)); GS=sparse(eye(mv*nv).*repmat(Svec,1,mv*nv)); 
% Inicialización de las matrices solución 
% Fuerza de los torbellinos de fronteras 
Gammab=zeros(mv*nv,ntimes); 
% Derivada en tiempo de las fuerzas de los torbellinos de fronteras 
Gammabdot=zeros(mv*nv,ntimes); 
% Fuerzas de los torbellinos de la estela 
Gammaw=zeros(mw*nv,1); 
Gammawt=zeros(mw*nv,ntimes); 
% Distribución de presiones 
p=zeros(mv*nv,ntimes); 
% Distribución de la sustentación 
l=zeros(mv*nv,ntimes); 
% Distribución del coeficiente de sustentación 
cl=zeros(mv*nv,ntimes); 
% Matriz de superficies para el calculo de p 
Svecmat=zeros(mv*nv,ntimes); 
% Inicio del bucle temporal: Calculo de las presiones y la sustentacion 
for it=1:ntimes 
    if it > 1 
        % Propagación de la estela 
        Gammaw=Pb*Gammab(:,it-1)+Pw*Gammaw; 
        Gammawt(:,it)=Gammaw; 
    end
    % Cálculo de la vorticidad de las fronteras 
    Gammab(:,it)=Ab\(-(Usalto)*sin(alpha)*nvec(:,3)-Aw*Gammaw); 
    % Cálculo de la derivada en el tiempo de la vorticidad de las fronteras 
    if it > 1 
        Gammabdot(:,it)=(Gammab(:,it)-Gammab(:,it-1))/dt; 
    else
        Gammabdot(:,1)=0; 
    end
    % Calculo de la distribucion de la sustentacion 
    l(:,it)=(rho*(Usalto)*cos(alpha)*Gy*Gammab(:,it)+ rho*GS*Gammabdot(:,it))*cos(alpha);
    % Calculo de la distribucion de presiones 
    Svecmat(:,it)=Svec(:); 
    p(:,it)=l(:,it)./(Svecmat(:,it)); 
end
% Sustentacion total 
L=sum(l,1); 
% Coeficiente de sustentacion 
Stotalala=sum(sum(S)); 
CL=L/(1/2*rho*(Usalto)^2*Stotalala); 
% Representacion del coeficiente de sustentacion total 
figure(1) 
plot(ti,CL,'color','black','LineWidth',2); 
xlabel('t(s)') 
ylabel('CL') 
title(['Coeficiente de sustentación total, CL(x,y,t), (mv = ', num2str(mv),', nv = ',num2str(nv), ')']) 
grid 
% Representacion de la última distribucion de la sustentacion 
figure(2) 
mesh(xp,yp,zp); 
hold on; 
% Sustentacion ultimo instante temporal 
lend=reshape(l(:,end),nv,mv)'; 
surf(xc,yc,lend); 
xlabel('Cuerda (metros)') 
ylabel('Envergadura (metros)') 
zlabel('l(x,y) (Newton)') 
title(['Distribución de la sustentación en t_f_i_n_a_l, l(x,y), (mv = ', num2str(mv),', nv = ',num2str(nv), ')']) 
grid on colorbar 
% Representacion de la última distribucion de presiones 
figure(3) 
mesh(xp,yp,zp); 
hold on 
% Presion ultimo instante temporal 
pend=reshape(p(:,end),nv,mv)'; 
h=surf(xc,yc,pend); 
xlabel('Cuerda (metros)') 
ylabel('Envergadura (metros)') 
zlabel('p(x,y) (Newton/m^2)') 
title(['Distribución de la presion en t_f_i_n_a_l, p(x,y), (nx = ', num2str(mv),', ny = ',num2str(nv), ')']) 
grid on 
colorbar

function [CF,CR,IF,IR,DF,DR]= EsquinasAla(LF,LR,b,cr) 
% CF: Coordenadas Centro - Frontal 
% CR: Coordenadas Centro - Posterior 
% IF: Coordenadas Izquierdo - Frontal 
% IR: Coordenadas Izquierdo - Posterior 
% DF: Coordenadas Derecho - Frontal 
% DR: Coordenadas Derecho - Posterior 
% Transformacion de angulos de flechas en grados a radianes 
LFrad = LF * pi/180; 
LRrad = LR * pi/180; 
% Cuerda en la punta de ala 
ct = cr - (b/2)*(tan(LFrad) - tan(LRrad)); 
% Coordenadas de las esquinas 
% "C": Central 
% "I": Izquierdo 
% "D": Derecho 
% "F": Frontal (ataque ala) 
% "R": Posterior (salida ala) 
% Coordenadas cuerda central 
CF = [ 0 , 0 ]; 
CR = [ cr , 0 ]; 
% Coordenadas cuerda extremo izquierdo(I) 
IF = [ (b/2)*tan(LFrad) , -b/2 ]; 
IR = [ (b/2)*tan(LFrad) + ct , -b/2 ]; 
% Coordenadas cuerda extremo derecho (D) 
DF = [ (b/2)*tan(LFrad) , +b/2 ]; 
DR = [ (b/2)*tan(LFrad) + ct , +b/2 ];
end

function [R] = CoordenadasPuntoInterior(xi,eta,LF,LR,b,cr) 
% R: Coordenadas R(x,y) del punto 
% xi: Coordenada 'x' adimensional con la cuerda 
% eta: Coordenada 'y' adimensional con la cuerda 
% LF: Ángulo de flecha frontal en grados (borde ataque) 
% LR: Ángulo de flecha posterior en grados (borde salida) 
% b: Envergadura en metros 
% cr: Cuerda en el encastre en metros 

% eta = 2y/b    % eta=0 (centro), 
                % eta=1 (punta semi-ala derecha)
                % eta=-1 (punta semi-ala izquierda) 
                % eta<0 (el punto se encuentra en el semi-ala izquierda) 
                % eta>0 (el punto se encuentra en el semi-ala derecha) 
%
% El punto en cuestión se encuentra a 100xi(%) del borde de ataque 
% para el valor correspondiente eta. 
% 
%               xi=0 (punto situado en el borde de ataque) 
%               xi=0.25 (punto situado en el c.a.) 
%               xi=1 (punto situado en el borde de salida)
%
% Obtenemos los puntos de las esquinas 
[CF,CR,IF,IR,DF,DR]=EsquinasAla(LF,LR,b,cr); 
if eta > 0 
    % CASO 1: PUNTO EN SEMI-ALA DERECHA, CUANDO eta > 0 
    R = (1 - eta) * xi * CR + eta * xi * DR + eta * (1 - xi) * DF + (1-eta)*(1-xi) * CF; 
elseif eta < 0 
    % CASO 2: PUNTO EN SEMI-ALA IZQUIERDA, CUANDO eta < 0 
    R = (1 + eta) * xi * CR + - eta * xi * IR + - eta * (1 - xi) * IF + (1+eta)*(1-xi) * CF; 
else
    % eta = 0 % CASO 3: PUNTO EN EL EJE x, CUANDO eta = 0 
    R = xi * CR + (1 - xi) * CF; 
end
end

function [fila,columna] = FilaColumnaNodo(k,py) 
% fila: Fila del nodo 
% columna: Columna del nodo 
% k: Numeracion del nodo 
% py: Número de nodos en direccion de la envergadura 
% Obtenemos los números naturales A y B tales que "k = A * py + B" 
A = floor(k/py); 
% Parte entera 
B = k - A * py; 
% Obtenemos la fila y la columna 
if B == 0 
    fila = A; 
    columna = py; 
else
    fila = A + 1; 
    columna = B; 
end
end

function [Mn] = MatrizNodos(nx,ny,LF,LR,b,cr) 
    % nx: Número de paneles en dirección x 
    % ny: Número de paneles en dirección y 
    % LF: Ángulo de flecha frontal en grados (borde ataque) 
    % LR: Ángulo de flecha posterior en grados (borde salida) 
    % b: Envergadura 
    % cr: Cuerda en la raíz 
    % La matriz Mn numera los nodos y asocia a cada uno sus coordenadas y otros 
    % datos 
    % Número de nodos en total
    % Puntos (nodos) en dirección x 
    px = nx + 1; 
    % Puntos (nodos) en dirección y 
    py = ny + 1; 
    % En total "P" puntos o nodos 
    P = px * py ; 
    % El nodo "1" es el situado en la esquina caracterizada por 
    % - Si el ala es recta y rectangular el que que tiene "y" más negativo 
    % - Si el ala es recta y rectangular el que que tiene "x" más positivo 
    % 
    % Con el eje "x" hacia abajo y el eje "y" hacia la derecha es el de la 
    % esquina inferior-izquierda. 
    % 
    % Los nodos se van numerando de izquierda a derecha y de abajo a arriba 
    % (y con ese criterio se numeran las columnas y filas (ver croquis) 
    % 
    % Ejemplo nx = 3, ny = 4, N=3*4 = 12 
    % px = 4, py = 5, P=4*5 = 20 
    % 
    % Matriz de nodos (numeración k=1,...20) 
    % Matriz de paneles [1]...[N=12]
    % 
    % f=4 16--------17---------18---------19----------20----------------(y) 
    %     |    [9]   |   [10]   |   [11]   |   [12]    | 
    % f=3 11--------12---------13---------14----------15 
    %     |    [5]   |    [6]   |    [7]   |    [8]    | 
    % f=2 6----------7----------8----------9----------10 
    %     |    [1]   |    [2]   |    [3]   |    [4]    | 
    % f=1 1----------2----------3----------4-----------5 
    %    c=1        c=2        c=3        c=4         c=5 
    %                           | 
    %                           | 
    %                           | 
    %                          (x) 
    % Configuración matriz P 
    % 
    % columna 1: Numeración del nodo (k) 
    % columna 2: fila del panel en la "matriz" de paneles (f) 
    % columna 3: columna del panel en la "matriz" de paneles (c) 
    % columna 4: Coordenada x del nodo k 
    % columna 5: Coordenada y del nodo k
    %
    % Inicialización matriz de paneles 
    Mn = zeros(P,5); 
    % Numeración de paneles 
    Mn(:,1) = 1:P; 
    for k=1:P 
        % Fila y columna 
        [f,c] = FilaColumnaNodo(k,py); 
        % Asignación a las columnas de P 
        Mn(k,2) = f; 
        Mn(k,3) = c; 
        % Coordenadas adimensionales del punto 
        xi = 1 - (f - 1)/nx; 
        eta = 2*(c - 1)/ny - 1; 
        % Cálculo de las coordenadas del punto 
        R = CoordenadasPuntoInterior(xi,eta,LF,LR,b,cr);
        x = R(1); 
        y = R(2); 
        Mn(k,4) = x; 
        Mn(k,5) = y; 
    end
end

function uvw=segmentotorbellino(x,y,z,x1,y1,z1,x2,y2,z2) 
% Cálculo de la velocidad inducida (u,v,w) en un punto (x,y,z) debido a un 
% segmento de torbellino de coordenadas (x1,y1,z1) a (x2,y2,z2) y fuerza 
% la unidad. 
% Distancia entre (x,y,z) y (x1,y1,z1) 
dist1=[x y z]-[x1 y1 z1]; 
% Distancia entre (x,y,z) y (x1,y1,z1) 
dist2=[x y z]-[x2 y2 z2]; 
% Cálculo del producto vectorial 
crossprod=cross(dist1,dist2); 
% Cálculo de la fracción 
fraction=crossprod/(crossprod*crossprod'); 
% Cálculo del producto escalar 
dotprod=dot(([x2 y2 z2]-[x1 y1 z1]),(dist1/sqrt(dist1*dist1')-dist2/sqrt(dist2*dist2'))); 
% Cálculo de la influencia total 
uvw=1/4/pi*fraction*dotprod; 
end
        