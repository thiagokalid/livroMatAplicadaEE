% matriz_admitancia realiza o c�lculo da matriz admit�ncia a partir da
% matriz imped�ncia das LTs e matriz suscpet�ncia das LTs.
%
%
%[matriz_Y] = matriz_admitancia(matriz_Z , matriz_YG)
%
% @param matriz_Z: matriz com valores da imped�ncia da LTs de tamanho
% N x N, onde N corresponde ao n�mero de barras.
% @param matriz_YG: matriz coluna com os valores das admit�ncias do 
% barramento para o terra. 
%
% @return matriz_Y: matriz quadrada (N x N) contendo os dados das
% admit�ncias das LTs conectando os barramentos (matriz_Y(n,m) = 0 significa que
% n�o h� LT entre os barramentos n e m).

function [matriz_Y] = matriz_admitancia(matriz_Z , matriz_YG)
% Todo o m�todo de c�lculo da matriz admit�ncia � detalhado melhor no livro
% do Ned Mohan Electric Power System: A first Course.
% 
% 
% -> Para elementos fora da diagonal principal, sua f�rmula resume-se ao
% inverso negativo da imped�ncia.
% 
% -> Para elementos dentro da diagonal pricipal, dever� ser somada a auto
% admit�ncia especificada matriz_YG e tamb�m todos elementos que tem
% conex�o (LT) com o barramento, ou seja, todos elementos de sua linha da matriz_Z.
% 
% * Caso n�o houver LT o valor da matriz_Z ser� zero, portanto n�o haver�
% contribui��o no c�lculo do somat�rio da admit�ncia.
% 

n = size(matriz_Z , 1);
matriz_Y = zeros(n,n);

for i = 1:n
    somatorio_linha = 0;
    for j = 1:n
        if(matriz_Z(i,j) ~= 0)
            % Primeiro ir� ser determinado todos os valores da linha que n�o
            % s�o da diagonal principal:
            matriz_Y(i,j) = -(1 / matriz_Z(i,j));
            somatorio_linha = somatorio_linha + (matriz_Y(i,j));
        end
    end
    % Ap�s determinado todos os valores, ou seja, chegado no final
    % da linha, ser� determinado o valor da diagonal principal:
    matriz_Y(i,i) = matriz_YG(i) - somatorio_linha;
end
end

