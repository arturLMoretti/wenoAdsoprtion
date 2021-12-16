function dCdt=RHS_exemplo_WENO(t,C,N,Deltaz,Pe,PSI1,PSI2,PSI3,PSI4,n)


%% Divis�o do vetor princiapal em X e Q
X1=C(1:N);
Q1=C(N+1:2*N);

% X2=C(2*N+1:3*N);
% Q2=C(3*N+1:4*N);

%% Defini��o das condi��es de contorno
CC1_1=Pe(1)*(X1(1)-1);
CC2_1=0;

%% C�lulas imagin�rias

% Como o WENO3 utiliza uma c�lula � direita e uma � esquerda, nas condi��es
% de contorno, uma aproxima��o para estes pontos n�o definidos deve ser
% realizada. Esta � uma recomenda��o das notas de aula do Shu, que mandei
% por e-mail. Outros m�todos tamb�m s�o v�lidos. 
X_antes_1=X1(1)-Deltaz*CC1_1;
X_depois_1=X1(N)+Deltaz*CC2_1;
X1_TEMP=[X_antes_1;X1;X_depois_1];

%% Obten��o do fluxo pelo WENO

% Mando o WENO3, que utiliza 3 pontos e o WENO5, que utliza 5 pontos. Os
% pontos adicionais do WENO5 est�o definidas na fun���
% 'weno_5_CC1_Neumann.m'. Para isso, utlizei extrapola�� linear, ou como j�
% li, o m�todo de aproxima��o da meia c�lula 

X1_z=weno_5_CC1_Neumann (X1,Deltaz,X_antes_1,X_depois_1);



%% Obten��o do fluxo por diferen�as finitas/ C�lculo da derivada segunda

% Pr�-alocando X_zz
X1_zz(1:N,1)=0;
for i=2:N+1
% Suprimindo a linha abaixo, o fluxo ser� dado pelo WENO, ao retirar o coment�rio no in�cio,
% o fluxo ser� aproximado por diferen�as finitas:

    %X_z(i-1)=1/(2*Deltaz)*(X_TEMP(i+1)-X_TEMP(i-1)); 
     
% Quando se considera a difus�o axial constante, pode-se mostrar que a
% discretiza��o em volumes finitos do fluxo difusivo pode ser aproximada
% pela derivada segunda

    X1_zz(i-1)=1/(Deltaz)^2*(X1_TEMP(i+1)-2*X1_TEMP(i)+X1_TEMP(i-1));
end

%% Defini��o da isoterma de equil�brio

Q1_EQ=isoterma_sips(PSI4(1),PSI3(1),n(1),X1);

%% Defni��o da quantidade m�dia adsorvida, pelo modelo LDF
dQ1dt=PSI2(1)*(Q1_EQ-Q1);

%% Defini��o do BM
dXdt=1/Pe(1)*X1_zz-X1_z-PSI1(1)*dQ1dt;

%% Unindo os vetors que fornecem o lado direito do sistema de EDOs obtidos da discretiza��o
dCdt=[dXdt;dQ1dt];
end