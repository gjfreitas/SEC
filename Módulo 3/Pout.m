clear all
close all
clc

% Parâmetros - Unidades SI
RHR = 1;        % 100%
xa = 0.22;      % 22 cm
xg = 0.18;      % 18 cm
gamma = 0.005;  %  m–1 
g0 = 0.15;      % m–1 
Psat = 1; %W

% Vamos variar 0% < Toc < 5% ( Nota: a expressão do Optimal Output Coupling apenas é valida para valores de Toc < 10 % 
Toc=linspace(0,0.05,100);
Pout=1/2*Psat*Toc.*((2*xg*g0)./(2*gamma*xa+Toc)-1);
Roc=1-Toc;

% Cálculo do valor ótimo do Toc do Optimal Output Coupling 
Toc_optimal = sqrt(4*xg*gamma*xa*g0)-2*gamma*xa;

% Solução exata de Rigrod (para low loss e low gain lasers)
R1 = RHR * (1-2*gamma*xa);  % lumped loss approximation para R1
Pout_e=(1-Roc).*Psat.*(g0*xg+log(sqrt(R1.*Roc)))./(1+sqrt(Roc/R1)-Roc-sqrt(R1*Roc));

figure(1)
hold on
plot(Roc,Pout,'b-','LineWidth',1,'MarkerSize',10)
plot(Roc,Pout_e,'r--','LineWidth',1,'MarkerSize',10)
xlabel('R_{oc}')
ylabel('P_{out} (W)')
legend('Optimal Output Coupling', 'Solução exata de Rigrod')
set(gca,'Xdir','reverse')



