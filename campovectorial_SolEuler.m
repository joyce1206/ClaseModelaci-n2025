%% Ejercicio dispersión de particulas con Euler y Runge Kutta 2
clear all
close all
clc
%% hasta aquí se grafica un campo vectorial
[x, y]= meshgrid(-2:.1:2); %Define el paso
%define X y Y como u y v
u=-y; v=x;
%u=y, v=sin(x);
%quiver(u,v), axis square
quiver(x,y,u,v)
%grid on
hold on
%% cálculo de Euler

delt=0.01;
xn=0.3;
yn=1.8;

plot(xn,yn,'om')
t=0;
hold on

for ii=1:500
    um=interp2(x,y,u,xn,yn)
    vm=interp2(x,y,v,xn,yn)

    k1x=delt*um;
    k1y=delt*vm;

    xn=xn+k1x;
    yn=yn+k1y;
    plot(xn,yn,'ob')
    pause(0.1)
end

%% Euler con varias particulas
%%%%%% Erupción instantánea %%%%%%

%Grafico del campo de viento
figure()
clf
quiver(x, y, u, v)
hold on

%Posición inicial de las particulas
XN = 1 + 0.01*randn(1,100); %Se modifica la posición aleatoriamente 
YN = 1 + 0.01*randn(1,100); %Se modifica la posición aleatoriamente 
plot(XN,YN,'m.')
%hold on

%Paso de tiempo
DelT = 0.01;

for ii=1:500
    %Velocidad en los puntos instantáneos
    UM = interp2(x,y,u,XN,YN);
    VM = interp2(x,y,v,XN,YN);
    
    %Calcula incremento con Euler
    k1x = UM;
    k1y = VM;
    
    %Agrega turbulencia a las particulas aleatoria
    npart = length(XN);
    XN = XN + DelT*k1x + 0.01*randn(1,npart);
    YN = YN + DelT*k1y + 0.01*randn(1,npart);

    clf
    quiver(x, y, u, v)
    hold on
    plot(XN,YN,'m.')
    axis([-2 2 -2 2])
    pause(0.1)
end
%% Pluma
figure()
clf
quiver(x, y, u, v)
hold on

XN = 1 + 0.01*randn(1,100);
YN = 1 + 0.01*randn(1,100);
DelT = 0.01;
plot(XN,YN,'m.')
%hold on

for ii=1:500
    %Velocidad en los puntos instantáneos
    UM = interp2(x,y,u,XN,YN);
    VM = interp2(x,y,v,XN,YN);
    
    %Calcula incremento con Euler
    k1x = DelT*UM;
    k1y = DelT*VM;
    
    %Agrega turbulencia a las particulas
    npart = length(XN);
    XN = XN + k1x + 0.01*randn(1,npart);
    YN = YN + k1y + 0.01*randn(1,npart);
    
    %Agrega más particulas durante el evento
    if ii < 20
        XN = [XN, 1.0 + 0.01*randn(1,10)];
        YN = [YN, 1.0 + 0.01*randn(1,10)];
    end
    
    %Grafica de posición iterativa
    clf
    quiver(x, y, u, v)
    hold on
    plot(XN,YN,'m.')
    axis([-2 2 -2 2])
    pause(0.1)
end
%% Con Runge Kutta 2 (a mano)
%k1 = h*f(x,y)
%k2 = h*f(x+h/2,y)

figure()
clf
quiver(x, y, u, v)
hold on

XN = 1 + 0.01*randn(1,100);
YN = 1 + 0.01*randn(1,100);
DelT = 0.01;
plot(XN,YN,'m.')
%hold on

for ii=1:100
    %Velocidad en los puntos instantáneos
    UM = interp2(x,y,u,XN,YN);
    VM = interp2(x,y,v,XN,YN);

    %Calcula incremento con Ruge Kutta 1
    k1x = UM;
    k1y = VM;

    %Calcula incremento con Ruge Kutta 2
    UM2 = interp2(x,y,u,XN+(3*k1x*DelT/4),YN+(3*k1x*DelT/4));
    VM2 = interp2(x,y,v,XN+(3*k1x*DelT/4),YN+(3*k1y*DelT/4));

    k2x = UM2;
    k2y = VM2;

    %Agrega turbulencia a las particulas
    npart = length(XN);
    XN = XN + (1/3)*(k1x + 2*k2x)*DelT + 0.01*randn(1,npart);
    YN = YN + (1/3)*(k1y + 2*k2y)*DelT + 0.01*randn(1,npart);
    
    %Agrega más particulas durante el evento
    if ii<20
        XN = [XN 1  + 0.001*randn(1,10)];
        YN = [YN 1  + 0.001*randn(1,10)];
    end
    
    %Grafica de posición iterativa
    clf
    quiver(x, y, u, v)
    hold on
    plot(XN,YN,'m.')
    pause(0.1)
end

%% Runge-Kutta2 con función advpartlb2
%Con una sola particula
%Posiciones iniciales
XN(1) = 1;
YN(1) = 0;

for k=1:100
    [XN YN] = advpartlb2(x,y,u,v,XN,YN,5,.1); % (Coord de la malla, particulas, pasos de tiempo)
    plot(XN, YN, '.r')
    axis([-2 2 -2 2])
    pause(0.01)
end

%% Con Runge Kutta 2 (hacia atras)
%k1 = h*f(x,y)
%k2 = h*f(x+h/2,y)

figure()
clf
quiver(x, y, u, v)
hold on

XN = 0.8 ;%+ 0.01*randn(1,100);
YN = 0.8 ;%+ 0.01*randn(1,100);
DelT = 0.01;
plot(XN,YN,'m.')
%hold on

for ii=1:100
    %Velocidad en los puntos instantáneos
    UM = interp2(x,y,u,XN,YN);
    VM = interp2(x,y,v,XN,YN);

    %Calcula incremento con Ruge Kutta 1
    k1x = UM;
    k1y = VM;

    %Calcula incremento con Ruge Kutta 2
    UM2 = interp2(x,y,u,XN+(k1x*DelT/2),YN+(k1y*DelT/2));
    VM2 = interp2(x,y,v,XN+(k1x*DelT/2),YN+(k1y*DelT/2));
    %UM2 = interp2(X1,Y1,U,XN+(3*DelT/4),YN+(3*k1x*DelT/4));
    %VM2 = interp2(X1,Y1,V,XN+(3*DelT/4),YN+(3*k1y*DelT/4));

    k2x = UM2;
    k2y = VM2;

    %Agrega turbulencia a las particulas
    npart = length(XN);
    XN = XN - (1/3)*(k1x + 2*k2x)*DelT ;%+ 0.01*randn(1,npart);
    YN = YN - (1/3)*(k1y + 2*k2y)*DelT ;%+ 0.01*randn(1,npart);
    
    %Agrega más particulas durante el evento
    if ii<20
        XN = [XN 0.8  + 0.001*randn(1,10)];
        YN = [YN 0.8  + 0.001*randn(1,10)];
    end
    
    %Grafica de posición iterativa
    clf
    quiver(x, y, u, v)
    hold on
    plot(XN,YN,'m.')
    pause(0.1)
end

%% Dibujando trayectoria con Euler
DelT = 0.01;

figure()
clf
quiver(x, y, u, v)
hold on

%Posicion inicial
XN = 0.3;
YN = 1.8;

plot(XN,YN,'om')
t=0;
hold on

tic; % Inicia el contador de tiempo
for ii=1:700
    UM = interp2(x,y,u,XN,YN);
    VM = interp2(x,y,v,XN,YN);
    
    k1x = UM;
    k1y = VM;
    
    XN = XN + DelT*k1x;
    YN = YN + DelT*k1y;
    plot(XN,YN,'.b')
    pause(0.01)
end
time_01 = toc; % Finaliza el contador de tiempo y guarda el valor
%% Euler paso de tiempo más pequeño
DelT = 0.001;

%Posicion inicial
XN = 0.3;
YN = 1.8;

plot(XN,YN,'om')
t=0;
hold on

tic; % Inicia el contador de tiempo
for ii=1:7000
    UM = interp2(x,y,u,XN,YN);
    VM = interp2(x,y,v,XN,YN);
    
    k1x = UM;
    k1y = VM;
    
    XN = XN + DelT*k1x;
    YN = YN + DelT*k1y;
    plot(XN,YN,'.r')
    pause(0.001)
end
time_001 = toc; % Finaliza el contador de tiempo y guarda el valor

%% Dibujando trayectoria con Runge-Kutta 2 (RK2)
DelT = 0.01;

% figure()
% clf
% quiver(x, y, u, v)
%hold on

% Posicion inicial
XN = 0.3;
YN = 1.8;

plot(XN, YN, 'om')
t = 0;
hold on

tic; % Inicia el contador de tiempo
for ii = 1:700
    % Velocidad en los puntos instantáneos
    UM = interp2(x, y, u, XN, YN);
    VM = interp2(x, y, v, XN, YN);

    % Calcula incremento con Runge-Kutta 1
    k1x = UM;
    k1y = VM;

    % Calcula incremento con Runge-Kutta 2
    UM2 = interp2(x, y, u, XN + (DelT / 2) * k1x, YN + (DelT / 2) * k1y);
    VM2 = interp2(x, y, v, XN + (DelT / 2) * k1x, YN + (DelT / 2) * k1y);

    k2x = UM2;
    k2y = VM2;

    % Actualiza las posiciones
    XN = XN + DelT * k2x;
    YN = YN + DelT * k2y;
    plot(XN, YN, '.g')
    pause(0.01)
end
time_rk = toc; 
%%
% Muestra los tiempos de ejecución
fprintf('Tiempo de ejecución con DelT = 0.01: %f segundos\n', time_01);
fprintf('Tiempo de ejecución con DelT = 0.001: %f segundos\n', time_001);
fprintf('Tiempo de ejecución runge Kutta 2 DelT = 0.01: %f segundos\n', time_rk);
