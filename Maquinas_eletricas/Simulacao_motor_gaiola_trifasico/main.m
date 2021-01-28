clear
clc

% Solu��o das simula��es do cap�tulo 1.
% Livro: "The Field Orientation Principle in Control of Induction Motors"
% Autor: Andrzej M. Trzynadlowski
%
% In�cio da se��o de entrada de dados do usu�rio:
%
% Dados da fonte trif�sica AC:
% V(t) = VLLrms * sin(2*pi*f*t + theta);
%
% Tens�o pico em volts:
VLLrms = 220;

% Frequ�ncia em hertz:
f = 60;

% Defasagem das fases:
theta = pi/2;

%
% Dados do motor:
%
% N�mero de polos do estator:
p = 6;

% Grau de escorregamento do motor:
s = 0.03;

% Resist�ncia do estator em ohm:
Rs = 0.288;

% Reat�ncia do estator em ohm:
Xls = 0.524;

% Resist�ncia do rotor em rela��o ao estator em ohm:
Rr = 0.158;

% Reat�ncia do rotor em rela��o ao estator em ohm:
Xlr = 0.279;

% Reat�ncia magnetizadora em ohm:
Xm = 15.457;

% Mom�nto de in�rcia do eixo do rotor em kilograma por metro quadrado:
Jm = 0.4;

% Mom�nto de in�cia da carga em kilograma por metro quadrado:
JL = Jm;


% Par�metros para solu��o num�rica das EDO do problema:
% Tempo inicial e final em segundos:
t0=0;
tf=5;
passo = 1e-6;

% O m�todo num�rico para solu��o das EDO ser� o de Euler.
% Fim da se��o de entrada de dados do usu�rio.

%%
% In�cio da se��o de c�lculos pr�vios antes do loop principal:
% Condi��es inicias das EDO a serem resolvidas:
% Velocidade �ngular mec�nica inicial em rad/s:
Wm_0 = 0;

w0 = p/2*Wm_0;

% Calculo da corrente inicial da equa��o diferencial:
[id_s0 , iq_s0 , id_r0 , iq_r0] = calc_idq0(VLLrms , Xlr , Xls, Xm , Rr , Rs , s);
i_dq0 = [id_s0 ; iq_s0 ; id_r0 ; iq_r0];

% Calculo do n�mero de itera��es da solu��o num�rica:
iteracoes = (tf-t0)/passo;

% C�lcilo da frequ�ncia angular ou velocidade �ngular s�ncrona em rad/s:
Wsync = 2 * pi * f; 

% Calculo das indut�ncias do motor em henry:
[Ls , Lm , Lr , L0] = calc_L(Xls , Xlr , Xm, Wsync);

% C�lculo do torque nominal do motor em newton metro:
T_rat = (p/2) * Lm * (iq_s0 * id_r0 - id_s0 * iq_r0);

% Vetores para armazenarem as vari�veis que ser�o plotadas:
Y = zeros(4,iteracoes);

% Definindo a matriz da transformada clark-parke:
T_abc2dq0 = sqrt(2/3) * ...
    [...
    cos(theta) cos(theta - 2/3*pi) cos(theta + 2/3*pi); ...
    -sin(theta) -sin(theta - 2/3*pi) -sin(theta + 2/3*pi); ...
    sqrt(2)/2 sqrt(2)/2 sqrt(2)/2 ...
    ];

% T_dq02abc = inv(T_abc2dq0)
T_dq02abc = sqrt(2/3) * ...
    [...
    cos(theta) -sin(theta) sqrt(2)/2 ;...
    cos(theta - 2/3*pi) -sin(theta - 2/3*pi) sqrt(2)/2 ;...
    cos(theta + 2/3*pi) -sin(theta + 2/3*pi) sqrt(2)/2 ;...
    ];

index = 1;

% Matriz que multiplica o vetor contendo as tens�es (dentro dos par�nteses) na equa��o 1.35:
matriz_A = [Lr 0 -Lm 0 ; 0 Lr 0 -Lm ;-Lm 0 Ls 0 ; 0 -Lm 0 Ls];

for t = t0:passo:tf
    % Fontes de tens�o CA:
    V_as = VLLrms*sqrt(2/3) * sin(Wsync*t + theta);
    V_bs = VLLrms*sqrt(2/3) * sin(Wsync*t + theta - 2/3*pi);
    V_cs = VLLrms*sqrt(2/3) * sin(Wsync*t + theta + 2/3*pi);
    
    % Vetor contendo as tr�s fases:
    V_abc_s = [V_as ; V_bs ; V_cs];

    
    V_dq0_s = T_abc2dq0 * V_abc_s;
    
    V_ds = V_dq0_s(1);
    V_qs = V_dq0_s(2);
    V_0s = V_dq0_s(3);
    
    % Como se trata de um motor tipo gaiola de esquilo, o rotor � curto-circuitado:
    V_dr = 0;
    V_qr = 0;
    
    V_dq = [V_ds ; V_qs ; V_dr  ; V_qr];


    % Calculando a derivada didt no ponto i_dq0 (equa��o 1.35):
    
    % Matriz que multiplica o vetor contendo as correntes (dentro dos par�nteses) na equa��o 1.35:
    matriz_B =[...
        -(Rs*Lr) w0*Lm^2 Rr*Lm w0*Lr*Lm ; ...
        -w0*Lm^2 -Rs*Lr -w0*Lr*Lm Rr*Lm ; ...
        Rs*Lm -w0*Ls*Lm -Rr*Ls -w0*Ls*Lr; ...
        w0*Ls*Lm Rs*Lm w0*Ls*Lr -Rr*Ls ...
        ];
    
    didt = 1/(L0^2) * (matriz_A * V_dq + matriz_B * i_dq0);
     
    % Solu��o num�rica da equa��o diferencial da corrente na forma matricial (equa��o 1.35):
    i_dq = i_dq0 + passo * didt;
    i_dq0 = i_dq;
    

    
    i_ds = i_dq(1);
    i_qs = i_dq(2);
    i_dr = i_dq(3);
    i_qr = i_dq(4);
    
    i_dq_s = [ i_ds ; i_qs ; 0]; 
    i_dq_r = [ i_dr ; i_qr ; 0];
    
    % C�lculo do torque (equa��o 1.49):
    
    T =  p/2 * Lm * (i_qs * i_dr - i_ds * i_qr);
    
    % Realizando a transformada dq0 para abc:
    
    i_abc_s = T_dq02abc * i_dq_s;
    i_abc_r = T_dq02abc * i_dq_r;
    
    
    % Definindo situa��o de carga constante durante a simula��o (simula�ao 1.1):
%     TL = T_rat*0.5;
    
    % Definindo situa��o de troca de carga durante simula��o (simula��o 1.2):
    if t < 3
        % Come�a com o torque da carga a 50 % do torque nominal da m�quina.
        TL = T_rat * 0.5;
    elseif t >= 3 && t < 4
        % Entre 3 e 4 segundos o torque da carga vai a 150 % do torque nominal da m�quina.
        TL = T_rat * 1.5;
    else
        % Depois dos 4 segundos o torque da carga vai a 100 % do torque nominal da m�quina.
        TL = T_rat * 1;
    end


    
    % Resolvendo a equa��o diferencial da velocidade angular mec�nica (equa��o 1.91):
    dWmdt = (T - TL)/(Jm + JL);
    Wm = Wm_0 + dWmdt * passo;
    Wm_0 = Wm;
    
    % Calculando a velocidade do rotor (equa��o 1.92):
    w0 = p/2*Wm;
    
    % Vetores que armazenam vari�veis que ir�o ser plotadas:
    Y(1,index) = Wm * 60/(2*pi);
    Y(2,index) = T;
    Y(3,index) = i_abc_s(1);
    Y(4,index) = i_abc_r(1);
    
    index = index + 1;
end

t = t0:passo:tf;


% Plotagem dos gr�ficos desejados:

subplot(2,2,1);
% figure(1);
plot(t,Y(1,:),"color","black");
grid on
ylabel("\omega_{M} in rpm");
xlabel("Time in s");
ylim auto
% xlim([2 5]);

subplot(2,2,2);
% figure(2);
plot(t,Y(2,:),"color","black");
grid on
ylabel("Torque in N \cdot m");
xlabel("Time in s");
ylim auto
% xlim([2 5]);

subplot(2,2,3);
% figure(3);
plot(t,Y(3,:),"color","black");
grid on
ylabel("i_{as} in A");
xlabel("Time in s");
ylim auto
% xlim([2 5]);

subplot(2,2,4);
% figure(4);
plot(t,Y(4,:),"color","black");
grid on
ylim auto
ylabel("i_{ar} in A");
xlabel("Time in s");
% xlim([2 5]);
