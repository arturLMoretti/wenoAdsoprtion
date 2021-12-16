function dCdt=RHS_exemplo_WENO(t,C,N,Deltaz,Pe,PSI1,PSI2,PSI3,PSI4,n)


%% Divisão do vetor princiapal em X e Q
X1=C(1:N);
Q1=C(N+1:2*N);

% X2=C(2*N+1:3*N);
% Q2=C(3*N+1:4*N);

%% Definição das condições de contorno
CC1_1=Pe(1)*(X1(1)-1);
CC2_1=0;

%% Células imaginárias

% Como o WENO3 utiliza uma célula à direita e uma à esquerda, nas condições
% de contorno, uma aproximação para estes pontos não definidos deve ser
% realizada. Esta é uma recomendação das notas de aula do Shu, que mandei
% por e-mail. Outros métodos também são válidos. 
X_antes_1=X1(1)-Deltaz*CC1_1;
X_depois_1=X1(N)+Deltaz*CC2_1;
X1_TEMP=[X_antes_1;X1;X_depois_1];

%% Obtenção do fluxo pelo WENO

% Mando o WENO3, que utiliza 3 pontos e o WENO5, que utliza 5 pontos. Os
% pontos adicionais do WENO5 estão definidas na funçãó
% 'weno_5_CC1_Neumann.m'. Para isso, utlizei extrapolaçã linear, ou como já
% li, o método de aproximação da meia célula 

X1_z=weno_5_CC1_Neumann (X1,Deltaz,X_antes_1,X_depois_1);



%% Obtenção do fluxo por diferenças finitas/ Cálculo da derivada segunda

% Pré-alocando X_zz
X1_zz(1:N,1)=0;
for i=2:N+1
% Suprimindo a linha abaixo, o fluxo será dado pelo WENO, ao retirar o comentário no início,
% o fluxo será aproximado por diferenças finitas:

    %X_z(i-1)=1/(2*Deltaz)*(X_TEMP(i+1)-X_TEMP(i-1)); 
     
% Quando se considera a difusão axial constante, pode-se mostrar que a
% discretização em volumes finitos do fluxo difusivo pode ser aproximada
% pela derivada segunda

    X1_zz(i-1)=1/(Deltaz)^2*(X1_TEMP(i+1)-2*X1_TEMP(i)+X1_TEMP(i-1));
end

%% Definição da isoterma de equilíbrio

Q1_EQ=isoterma_sips(PSI4(1),PSI3(1),n(1),X1);

%% Defnição da quantidade média adsorvida, pelo modelo LDF
dQ1dt=PSI2(1)*(Q1_EQ-Q1);

%% Definição do BM
dXdt=1/Pe(1)*X1_zz-X1_z-PSI1(1)*dQ1dt;

%% Unindo os vetors que fornecem o lado direito do sistema de EDOs obtidos da discretização
dCdt=[dXdt;dQ1dt];
end