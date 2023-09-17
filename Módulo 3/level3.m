clear all
close all
clc

y0=[1 0 0]; % Condição inicial para a população, em t=0 toda se encontra no nível 1
TSPAN = [0 0.05]; % Periodo de tempo para a simulação
[T,Y] = ode15s(@rate_eq,TSPAN,y0);

figure(1) % Niveis 1 e 3
hold on
plot(T,Y(:,1),'b-','LineWidth',1,'MarkerSize',10)
plot(T,Y(:,3),'r-','LineWidth',1,'MarkerSize',10)
legend('Nível 1', 'Nível 3')
xlabel('t')
ylabel('População')

figure(2) % Todos os Niveis
hold on
plot(T,Y(:,1),'b-','LineWidth',1,'MarkerSize',10)
plot(T,Y(:,2),'k-','LineWidth',1,'MarkerSize',10)
plot(T,Y(:,3),'r-','LineWidth',1,'MarkerSize',10)
legend('Nível 1', 'Nível 2', 'Nível 3')
axis([0 0.05 -0.1 1.1])
xlabel('t')
ylabel('População')

% Verificar que N(1)+N(2)+N(3) é sempre constante e igual a 1
somatorio = Y(:,1)+Y(:,2)+Y(:,3);

function dy = rate_eq(t,y)
    sigma = 10E-17; % cm^2 - Absorption cross-section
    Pd = 20; % W/cm^2 - Power density
    lambda = 600; % nm - Excitation wavelength
    h = 6.62 * 10E-27; % egr.s - Planck Constant
    W23 = 10E7; % s^-1 - 2-> 3 rate
    c = 29979245800; % cm/s - speed of light

    % Varia Tau entre (1.5,0.1,0.17)E-3 para obter todos os gráficos apresentados neste projeto
    Tau = 1.5 * 10E-3; % s - laser decay lifetime

    phi = ( sigma*Pd*lambda)/(h*c); % s^-1 - pumping rate

    dy = zeros(3,1);
    dy(1) = -(phi*y(1)) + ((1 / Tau) * y(3));
    dy(2) = (phi * y(1)) - (W23 * y(2));
    dy(3) = (W23 * y(2)) - ((1 / Tau) * y(3));
end
