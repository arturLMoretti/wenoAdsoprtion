function fluxo=weno_3_CC1_Neumann (u,dx,u_antes,u_depois)
% Retorna o fluxo calculado pelo WENO3. Considerando condi��es de Neumann
% nas duas extremidades do intervalo de resolu��o. 
% Entrada:
%   u -> vetor cujo fluxo deve ser calculado, cont�m N elementos
%   dx -> Tamanho da c�lula, cont�m 1 elemento
%   F1 -> Lado direito da condi��o de contorno da extremidade x=0, quando
%   escrita da forma: 
%       du/dx|_{x=0}= F1 (u,x)
%       Cont�m 1 elemento
%   F2 -> Lado direito da condi��o de contorno da extremidade x=1, quando
%   escrita da forma: 
%       du/dx|_{x=1}= F2 (u,x)
%       Cont�m 1 elemento

u= reshape(u,1,[]);
u=[u_antes,u,u_depois];

f=u;
dfdu=ones(size(f));

a=max(abs(dfdu)); f_mais=0.5*(f+a*u); f_menos=circshift(0.5*(f-a*u),[0,-1]);

%% Right Flux
% Choose the positive fluxes, 'v', to compute the left cell boundary flux:
% $u_{i+1/2}^{-}$
vm  = circshift(f_mais,[0 1]);
vp  = circshift(f_mais,[0 -1]);

% Polynomials
p0n = (-vm + 3*f_mais)/2;
p1n = ( f_mais  + vp )/2;

% Smooth Indicators, Beta factors
B0n = (vm-f_mais).^2; 
B1n = (f_mais-vp).^2;

% Constants
d0n = 1/3; d1n = 2/3; epsilon = 1E-6;

% Alpha weights 
alpha0n = d0n./(epsilon + B0n).^2;
alpha1n = d1n./(epsilon + B1n).^2;
alphasumn = alpha0n + alpha1n;

% ENO stencils weigths
w0n = alpha0n./alphasumn;
w1n = alpha1n./alphasumn;

% Numerical Flux at cell boundary, $u_{i+1/2}^{-}$;
hn = w0n.*p0n + w1n.*p1n;

%% Left Flux 
% Choose the negative fluxes, 'u', to compute the left cell boundary flux:
% $u_{i-1/2}^{+}$ 
um  = circshift(f_menos,[0 1]);
up  = circshift(f_menos,[0 -1]);

% Polynomials
p0p = ( um + f_menos )/2;
p1p = (3*f_menos - up)/2;

% Smooth Indicators, Beta factors
B0p = (um-f_menos).^2; 
B1p = (f_menos-up).^2;

% Constants
d0p = 2/3; d1p = 1/3; epsilon = 1E-6;

% Alpha weights 
alpha0p = d0p./(epsilon + B0p).^2;
alpha1p = d1p./(epsilon + B1p).^2;
alphasump = alpha0p + alpha1p;

% ENO stencils weigths
w0p = alpha0p./alphasump;
w1p = alpha1p./alphasump;

% Numerical Flux at cell boundary, $u_{i-1/2}^{+}$;
hp = w0p.*p0p + w1p.*p1p;

fluxo= (hp-circshift(hp,[0 1])+hn-circshift(hn,[0 1]))'/dx;

fluxo(1)=[];
fluxo(end)=[];

end