% Channel gain LOS
clear all
close all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Room Specs
Xmax = 6;               % Max width
Ymax = 6;               % Max length
h = 3;                  % Ceiling height
dx = 0.06;              % spatial resolution for ground dy=dx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TX and RX Specs
hpa = 60;               % TX half power angle
m = - log(2)/log(cosd(hpa));
FOV = 90;
R = 0.7;                % RX responsivity
Apd = 100e-6;           % area of the photo-detector
Pt = 1;                 % transmitted power
% TX positions
[XT, YT] = meshgrid([-2 2],[-2 2]);
NT = sqrt(prod(size(XT)));
nt = [0 0 -1];
% RX positions
[X, Y] = meshgrid(-Xmax/2+dx/2:dx:Xmax/2-dx/2,-Ymax/2+dx/2:dx:Ymax/2-dx/2);
NR = sqrt(prod(size(X)));
nr = [0 0 1];
% LOS contribution
Pr_LOS_TX = cell(NT,NT);
CLOS = (m+1)*R*Apd*Pt/(2*pi);
for t1 = 1:NT
    for t2 = 1:NT
        Pr = zeros(size(X));
        Tx_pos = [XT(t1,t2) YT(t1,t2) h];
        for r1 = 1:NR
            for r2 = 1:NR
                Rx_pos = [X(r1,r2) Y(r1,r2) 0];
                d = Rx_pos - Tx_pos;
                dn = sqrt(d*d.');
                cos_phi = (nt*d.')/dn;
                cos_psi = -(nr*d.')/dn;
                if acos(cos_psi) < FOV
                    Pr(r1,r2) = CLOS*((cos_phi)^m)*cos_psi/(dn^2);
                end
            end
        end
        Pr_LOS_TX{t1,t2} = Pr;
    end
end
figure(1)
Pr_LOS = PlotChannel(X,Y,NT,Pr_LOS_TX,'P_r_L_O_S (W)');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Auxiliary functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function M = PlotChannel(X,Y,NT,C,str)
% Makes a surface plot of the channel DC response for all transmitters
% X, Y specify the receiving plane grid
% NT is the square root of the number of transmitters in the room. Number
% of transmitters has to be a square number.
% C is a cell array, with the individual response due to each transmitter
% str is a string specifwying the type of response being ploted
    M = zeros(size(X));
    for t1 = 1:NT
        for t2 = 1:NT
            M = M + C{t1,t2};
        end
    end
    surfc(X,Y,M)
    shading interp
    xlabel('x (m)')
    ylabel('y (m)')
    zlabel(str)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%