%Solução da EDO:                 d2vc    1    dvc   1
%                           C  * ----- + - * ---- + -  * vc = 0
%                                 dt2    R    dt    L
%
% Condições inicias:  vc(0) = 5, vc'(0) = 0;
% Constantes: C = 1 F ; L = 1 H ; R = 1 Ohm;
%
%%

clear all
close all
clc


syms vc(t) vc1 R L C V0 t s


%%
%Este é o bloco das variáveis definidas pelo usuário:

c =1; % Valor do capacitor em Farad;
r =1; % Valor do resistor em Ohm;
l =1; % Valor do indutor em Henry;

%Intervalo de tempo da amostragem:
t0 = 0; % Tempo inicial;
tf = 15; % Tempo final;

% Definindo o tipo da fonte de tensão:
w = 377;
V(t) = 1200 * sin(w*t);

% Definindo as condições inicias:
Vc = 5; % Carga inicial no capacitor.
DVcDt = 0; % Derivada da tensão no capacitor.

%Tamanho do degrau das soluções numéricas:
step = 0.005;




%%
%Solução simbólica usando transformada de Laplace:

tic

% Definindo as derivadas simbólicas de ordem 1 e 2;
dvcdt = diff(vc,t);
d2vcdt2 = diff(dvcdt,t);

% Definindo as constantes simbolicamente e seus respectivos valores numéricos:
variavel_simbolica_da_constantes = [R L C];
valor_numerico_das_constantes =    [r l c];

% Definindo as condições iniciais simbolicamente e seus respectivos valores numéricos
variavel_simbolica_das_condicao_iniciais = [vc(0) dvcdt(0) d2vcdt2(0)];
valor_numerico_da_condicoes_inicial = [Vc DVcDt 0];

% Fazendo a transformada de laplace:
equation_t = V(t) == d2vcdt2 * L + dvcdt* R + vc * (1/C); % Equação diferencial em t;
equation_s = laplace(equation_t, t, s); % Equação algébrica em s;

% Simplificando e resolvendo a equação, agora algébrica, em s:
equation_s_sub_1 =  subs(equation_s ,[laplace(vc(t), t, s)],[vc1]);
equation_s_sub_2 = subs(equation_s_sub_1,variavel_simbolica_da_constantes,valor_numerico_das_constantes); % Substituindo os valores das constantes;
equation_s_sub_3 = subs(equation_s_sub_2, variavel_simbolica_das_condicao_iniciais, valor_numerico_da_condicoes_inicial); % Substituindo os valores das condições inicias;
solution_s = solve(equation_s_sub_3,vc1); % Resolvendo a equação em s, em função de vc1;

% Fazendo a transformada de laplace inversa para obter a solução em t da equação diferencial proposta:
solution_t = ilaplace(solution_s,s,t); % Transformando de s para t a solução;
simplified_solution_t = simplify(solution_t); % Simplificando a solução;

disp("Equação diferencial: ");
pretty(equation_t);

disp("Solução da equação diferencial simbolicamente: ");
pretty(simplified_solution_t);

disp(sprintf('Tempo de execução do método Simbolico: %f s',toc));


subplot(2,2,1);
fplot(simplified_solution_t,[t0 tf]);
grid on;
title("Solução simbólica ordem 2: Tensão no capacitor Vc(t)");
xlabel("t(s)");
ylabel("Vc(V)");


%%
% Formando a função numérica a partir da simbólica para aplicação nos métodos numéricos:

%Manipulações algébricas para isolar a derivada segunda:
segunda_derivada = isolate(equation_t , d2vcdt2);
simbolica_d2vcdt2_1 = rhs(segunda_derivada);
simbolica_d2vcdt2_2 = subs(simbolica_d2vcdt2_1,variavel_simbolica_da_constantes,valor_numerico_das_constantes);

valores = [t vc dvcdt];
syms dydt d2ydt2
valores_novos = [t dydt d2ydt2];

simbolica_d2vcdt2_3 = subs(simbolica_d2vcdt2_2,valores,valores_novos);
funcao_numerica = matlabFunction(simbolica_d2vcdt2_3,'Vars',[d2ydt2 t dydt ]);


%%
% Solução numérica usando o método de Euler:

tic
clear dvcdt d2vcdt2 % Reiniciando as derivadas pois agora iremos usar a função numérica obtida;

n = (tf - t0) / step; % Calculando número de iterações necessárias;
t = 0; % Variável tempo para plotagem;

v1_i = valor_numerico_da_condicoes_inicial(1);
v2_i = valor_numerico_da_condicoes_inicial(2);

dvcdt =@(vc)(vc); % Primeira derivada;
d2vcdt2 =funcao_numerica; % Segunda derivada;

subplot(2,2,2);

hold on
grid on

for i = 1:n
    v1_f = v1_i + step * dvcdt(v2_i);
    v1_i = v1_f;

    v2_f = v2_i + step * d2vcdt2(dvcdt(v2_i),t, v1_i); % Tensão no indutor divido pela auto-indutâcia
    v2_i = v2_f;
    
    plot(t, v1_i,'.', 'color', 'blue');
    t = t + step;
end


title("Euler ordem 2: Tensão no capacitor Vc(t)");
xlabel("t(s)");
ylabel("Vc(v)");

disp(sprintf('Tempo de execução do método Euler: %f s',toc));
%%
%Solução numérica usando o método trapezoidal:

tic
t = 0;
v1_i = valor_numerico_da_condicoes_inicial(1);
v2_i = valor_numerico_da_condicoes_inicial(2);

subplot(2,2,3);

hold on;
grid on;

for i=1:n

    v1 = v1_i + step * dvcdt(v2_i);
    v2 = v2_i + step * d2vcdt2( dvcdt(v2_i),t,v1_i);

    t=t+step;
    loopcount_1=0;
    diff=1;
    
    while abs(diff)>.5
        loopcount_1 = loopcount_1 +1;
        v1_final = v1_i + step *(dvcdt(v2_i)+dvcdt(v2_i))/2;
        diff= v1 - v1_final;
        v1 = v1_final;
    end
    
    loopcount_2=0;
    diff=1;
    
    while abs(diff)> .5

        loopcount_2 = loopcount_2 + 1;
        v2_final = v2_i + step *(d2vcdt2( dvcdt(v2_i),t,v1_i)+ d2vcdt2(dvcdt(v2),t,v1))/2;
        diff= v2 - v2_final;
        v2 = v2_final;
        
    end
    


plot(t,v1,'.','color','red');

v1_i = v1;
v2_i = v2;
end

title('Método Trapezoidal: Tensão no capacitor Vc(t)');
xlabel('t  (s)');
ylabel('Vc  (V)');


disp(sprintf('Tempo de execução do método Trapezoidal: %f s',toc));

%%
%Solução numérica usando o método Runge Kutta 2:

tic
t = 0;
v1_i = valor_numerico_da_condicoes_inicial(1);
v2_i = valor_numerico_da_condicoes_inicial(2);

subplot(2,2,4);

hold on
grid on
for i = 1:n
    k1_v1 = dvcdt(v2_i);
    k2_v1 = dvcdt((step * k1_v1) + v2_i);
    
    v1_f = v1_i + (step/2) * (k1_v1 + k2_v1);
    
    k1_v2 = d2vcdt2( dvcdt(v2_i),t,v1_i);
    k2_v2 = d2vcdt2( dvcdt(v2_i) + (step * k1_v2),t, v1_i + (step * k1_v2));
    
    v2_f = v2_i + (step/2) * (k1_v2 + k2_v2);
    
    v1_i = v1_f;
    v2_i = v2_f;
    
    plot(t, v1_i,'.', 'color', 'blue');
    t = t + step;
end

title("Runge Kutta ordem 2: Tensão no capacitor Vc(t)");
xlabel("t(s)");
ylabel("Vc(v)");

disp(sprintf('Tempo de execução do método Runge Kutta 2: %f s',toc));

ax(1) = subplot(2,2,1);
ax(2) = subplot(2,2,2);
ax(3) = subplot(2,2,3);
ax(4) = subplot(2,2,4);
linkaxes(ax,'x');
linkaxes(ax,'y');
