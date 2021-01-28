% contBarramentos conta quantos barramentos tipo PV e PQ existem no
% sistema.
%
%
% [nPV , nPQ] = contBarramentos(Tipo_barra)
%
%@param Tipo_barra: Vetor linha contendo as especifica��es de qual tipo �
%cada barramento.
%
% @param nPV: n�mero escalar inteiro que cont�m o quantidade de barramentos PV;
% @param nPQ: n�mero escalar inteiro que cont�m o quantidade de barramentos PQ;


function [nPV , nPQ] = contBarramentos(Tipo_barra)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
num_barramentos = size(Tipo_barra, 2); % Como � matriz coluna, o n�mero de linhas corresponde ao numer ode elementso.
nPV = 0;
nPQ = 0;

for i = 1:1:num_barramentos
    switch( Tipo_barra(i))
        case 1 % Slack
            
        case 2 % PV
            nPV = nPV + 1;
        case 3 % PQ
            nPQ = nPQ + 1;
    end
end
end

