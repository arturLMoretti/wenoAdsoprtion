clc
clear all



set(0,'DefaultAxesFontSize', 14,'DefaultAxesFontWeight','bold','DefaultAxesLineWidth',1.5);
set(0,'DefaultLegendFontSize',14,'DefaultLegendLocation','SouthEast','DefaultLegendFontSizeMode','manual');
set(0,'DefaultLineLineWidth',1.5);
set(0,'DefaultAxesGridLineStyle','-');


N=200;
Deltaz=1/N;
z=0:Deltaz:1;

% O problema se torna cada vez mais stiff quanto maior o Pe. Neste casos,
% as diferen�as finitas ainda d�o conta, mas muitas oscila��es s�o
% obsrvadas. Quando Pe=1000, j� � poss�vel observar isso. Quanto maior Pe,
% nem o WENO d� conta, mas acredito que da� as CCs devem ser modificadas

Pe_1=100;
PSI1_1=0;
PSI2_1=0;
PSI3_1=1;
PSI4_1=1;
n=[1e-1;1; 10];
NN=length(n);
t=linspace(0,10,5000);
C0=zeros(2*N,1);

for i=1:NN
    LEGENDA(i,1)={strcat('n','{ }','=','{ }',num2str(n(i)))};
end

figure
for i=1:NN
[t,C]=ode15s(@(t,C)RHS_exemplo_WENO(t,C,N,Deltaz,Pe_1,PSI1_1,PSI2_1,PSI3_1,PSI4_1,n(i)),t,C0);
plot(t,C(:,N))
xlabel 't^*'
ylabel 'X_i'
hold on
end
title 'Efeito do Par�metro n nas curvas de ruptura'
legend(LEGENDA)