% rbc_cai

% Miguel Alvarado

clear global;
clear all;
clc;
% cd D:\matlab_my_files\malvarado_programs\ejercicio2_b;
format compact;
format short;

global p_beta p_delta p_alfa p_sigma p_v p_chi p_rho ...
       c_ss k_ss z_ss n_ss y_ss i_ss w_ss R_ss q_ss lambda_ss;

% parametros

p_beta  = 0.9825;
p_delta = 0.025;
p_alfa  = 1/3;
p_sigma = 2;
p_v     = 1;
p_chi   = 1;
p_rho   = 0.9;

% Calculo Estado EStacionario (ver archivo ss_for_seq_rbc_cai.m)

run ss_for_seq_rbc_cai.m

% Para el sistema log-linealizado colapsado (ver archivo pdf, Sistema 2), 
% denotamos los siguientes par�metros

a_1 = -p_sigma/(p_v+p_alfa); %ok
a_2 = p_alfa/(p_v+p_alfa); %ok
a_3 = 1/(p_v+p_alfa); %ok

b_1 = y_ss/i_ss; %ok
b_2 = c_ss/i_ss; %ok

d_1 = ((1-p_alfa)*b_1*a_1)-b_2; %ok
d_2 = (p_alfa*b_1)+((1-p_alfa)*b_1*a_2); %ok
d_3 = b_1+((1-p_alfa)*b_1*a_3); %ok

% f_1 = -(1+p_beta*p_chi)/p_delta;
% f_2 = ((1+p_beta*p_chi)*(1-p_delta))/p_delta;

f_3 = p_beta*p_chi*d_1; %ok
f_4 = p_beta*p_chi*d_2; %ok
f_5 = p_beta*p_chi*d_3; %ok

g_1 = (1-p_alfa)*a_1; %ok
g_2 = (p_alfa-1)+((1-p_alfa)*a_2); %ok
g_3 = 1+((1-p_alfa)*a_3); %ok

h_1 = -(p_sigma-(p_beta*R_ss*g_1)); %ok
h_2 = p_beta*R_ss*g_2; %ok
h_3 = p_beta*R_ss*g_3; %ok
h_4 = p_beta*(1-p_delta); %ok

e_1 = y_ss/k_ss; %ok
e_2 = c_ss/k_ss; %ok

r_1 = ((1-p_alfa)*e_1*a_1)-e_2; %ok
r_2 = (1-p_delta)+(p_alfa*e_1)+((1-p_alfa)*e_1*a_2); %ok
r_3 = e_1+((1-p_alfa)*e_1*a_3); %ok

p_1 = -(1+p_beta)*p_chi*d_1; %ok
p_2 = -(1+p_beta)*p_chi*d_2; %ok
p_3 = -(1+p_beta)*p_chi*d_3; %ok


% modelo matricial;

A1 = [-p_sigma 1 0 0 0;
      p_1 1 p_chi p_2 p_3;
      d_1 0 0 d_2 d_3;
      r_1 0 0 r_2 r_3;
      0 0 0 0 p_rho];
  
A2 = [h_1 h_4 0 h_2 h_3;
      -f_3 0 0 -f_4 -f_5;
      0 0 1 0 0;
      0 0 0 1 0;
      0 0 0 0 1];

% formal estructural  
 
A1_inv = inv(A1);
A = A1_inv*A2;

%B = A1_inv*A3;

[Q,F] = jordan(A); %F matriz de autovalores de A, Q matriz de autovectores de A
[U,T] = schur(A);  %T matriz de autovalores de A, U matriz de autovectores de A

variables = 5;
vcontrol = 2; %variablles libres {c_(t), q_(t)}, I es forward looking

root_g_1=0; % contador
root_l_1=0; % contador

for i=1:variables
if abs(F(i,i))>=1
    root_g_1=root_g_1 + 1;
else
    root_l_1=root_l_1 + 1;
end
end;

if root_l_1==vcontrol
    display('SISTEMA ESTABLE'), F
else
    display('SISTEMA INESTABLE'), break
end

% Funciones de Politica (FP) & IRF

QQ=Q^-1;

q=zeros(root_l_1,variables); % vector donde recoger� los coeficientes para armar la FP-IRF

j=1; % contador

for i=1:variables
    if abs(F(i,i))<1
        q(j,:)=QQ(i,:);
        j=j+1;
    end
end

% Matriz de autovectores de raiz estable

q;

% vectores donde recoger� la trayectoria de las variables del modelo

z_tray = zeros(50,1);
c_tray = zeros(50,1);
k_tray = zeros(51,1);
I_tray = zeros(50,1);
q_tray = zeros(50,1);

% Din�mica del shock de 1% en el AR(1) del modelo
% de q, se define los siguiente par�metros:

om_1 = q(1,1)-((q(1,2)*q(2,1))/q(2,2)); %omega en la ecuacion de transici�n de c
om_2 = ((q(1,2)*q(2,3))/q(2,2))-q(1,3);
om_3 = ((q(1,2)*q(2,4))/q(2,2))-q(1,4);
om_4 = ((q(1,2)*q(2,5))/q(2,2))-q(1,5);

mu_1 = -(q(2,1)/q(2,2)); %mu en la ecuacion de transici�n de q
mu_2 = -(q(2,3)/q(2,2));
mu_3 = -(q(2,4)/q(2,2));
mu_4 = -(q(2,5)/q(2,2));

z_tray(1,1) = 1; % Shock en t = 1
k_tray(1,1) = 0; % Cero al ser variable predeterminada (estado) en t = 1.
I_tray(1,1) = 0; % Cero, dado que I_{t} = i_{t-1}, en t = 1, i_{0} no esta desviado de su estado estacionario.
c_tray(1,1) = ((om_2/om_1)*I_tray(1,1)) + ((om_3/om_1)*k_tray(1,1)) + ((om_4/om_1)*z_tray(1,1)); % Ecuaci�n 21
q_tray(1,1) = (mu_1*c_tray(1,1)) + (mu_2*I_tray(1,1)) + (mu_3*k_tray(1,1)) + (mu_4*z_tray(1,1)); % Ecuaci�n 22

% Del sistema final

for i=2:50
    z_tray(i,1) = p_rho*z_tray(i-1,1); % Ecuaci�n 25 
    k_tray(i,1) = (r_1*c_tray(i-1,1)) + (r_2*k_tray(i-1,1)) + (r_3*z_tray(i-1,1)); % Ecuaci�n 23
    I_tray(i,1) = (d_1*c_tray(i-1,1)) + (d_2*k_tray(i-1,1)) + (d_3*z_tray(i-1,1)); % Ecuaci�n 24
    c_tray(i,1) = ((om_2/om_1)*I_tray(i,1)) + ((om_3/om_1)*k_tray(i,1)) + ((om_4/om_1)*z_tray(i,1)); % Ecuaci�n 21
    q_tray(i,1) = (mu_1*c_tray(i,1)) + (mu_2*I_tray(i,1)) + (mu_3*k_tray(i,1)) + (mu_4*z_tray(i,1)); % Ecuaci�n 22
end

% calculamos un periodo adicional para el stock de capital y nuestra variable auxiliar.

k_tray(51,1) = (r_1*c_tray(50,1)) + (r_2*k_tray(50,1)) + (r_3*z_tray(50,1));
I_tray(51,1) = (d_1*c_tray(50,1)) + (d_2*k_tray(50,1)) + (d_3*z_tray(50,1));

% calculamos el resto de las trayectorias:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% vectores donde recoger� la trayectoria de las variables restantes del modelo
% y recuperar� la inversi�n desde mi variable auxiliar

lambda_tray = zeros(50,1);
n_tray  = zeros(50,1);
y_tray  = zeros(50,1);
w_tray = zeros(50,1);
R_tray = zeros(50,1);
i_tray  = zeros(50,1);

for i=1:50
    lambda_tray(i,1) = (-1)*p_sigma*c_tray(i,1); % Ecuacion 1
    n_tray(i,1) = (a_1*c_tray(i,1)) + (a_2*k_tray(i,1)) + (a_3*z_tray(i,1)); % Ecuacion 1 y 5 en 2
    w_tray(i,1) = z_tray(i,1) + p_alfa*k_tray(i,1) - p_alfa*n_tray(i,1); % Ecuacion 5
    R_tray(i,1) = z_tray(i,1) + (p_alfa - 1)*k_tray(i,1) + (1 - p_alfa)*n_tray(i,1); % Ecuacion 6
    y_tray(i,1) = z_tray(i,1) + p_alfa*k_tray(i,1) + (1 - p_alfa)*n_tray(i,1); % Ecuacion 7
    i_tray(i,1) = I_tray(i+1,1); % Ecuacion 11
end


% Ajusto el lag que realiza dynare en el stock de capital
k_tray = [k_tray(2:51,1)];

% Variables para los gr�ficos.
time = [0:1:49];
l_cero = zeros(50,1);

% Grafico de las IRF's

% Trayectoria tecnolog�a: z
subplot(3,4,1), plot(time,z_tray,'b-',time,l_cero,'r-','Linewidth',1.0), axis([0 49 -0.05 max(z_tray)+0.05]), title('Tecnolog�a: z')

% Trayectoria Consumo: c
subplot(3,4,2), plot(time,c_tray,'b-',time,l_cero,'r-','Linewidth',1.0), axis([0 49 -0.05 max(c_tray)+0.05]), title('Consumo: c')

% Trayectoria Stock de Capital: k
subplot(3,4,3), plot(time,k_tray,'b-',time,l_cero,'r-','Linewidth',1.0), axis([0 49 -0.05 max(k_tray)+0.05]), title('Stock de Capital: k')

% Trayectoria Trabajo: n
subplot(3,4,4), plot(time,n_tray,'b-',time,l_cero,'r-','Linewidth',1.0), axis([0 49 min(n_tray)-0.05 max(n_tray)+0.05]), title('Trabajo: n')

% Trayectoria Inversi�n: i
subplot(3,4,5), plot(time,i_tray,'b-',time,l_cero,'r-','Linewidth',1.0), axis([0 49 min(i_tray)-0.05 max(i_tray)+0.05]), title('Inversi�n: i')

% Trayectoria Producto: y
subplot(3,4,6), plot(time,y_tray,'b-',time,l_cero,'r-','Linewidth',1.0), axis([0 49 min(y_tray)-0.05 max(y_tray)+0.05]), title('Producto: y')

% Trayectoria Salarios: w
subplot(3,4,7), plot(time,w_tray,'b-',time,l_cero,'r-','Linewidth',1.0), axis([0 49 -0.05 max(w_tray)+0.05]), title('Salarios: w')

% Trayectoria Costo del Capital: R
subplot(3,4,8), plot(time,R_tray,'b-',time,l_cero,'r-','Linewidth',1.0), axis([0 49 min(R_tray)-0.05 max(R_tray)+0.05]), title('Costo del Capital: R')

% Trayectoria Lambda: Lambda
subplot(3,4,9), plot(time,lambda_tray,'b-',time,l_cero,'r-','Linewidth',1.0), axis([0 49 min(lambda_tray)-0.1 0.05]), title('Lambda: Lambda')

% Trayectoria q de Tobin: q
subplot(3,4,10), plot(time,q_tray,'b-',time,l_cero,'r-','Linewidth',1.0), axis([0 49 min(q_tray)-0.05 max(q_tray)+0.05]), title('q de Tobin: q')

