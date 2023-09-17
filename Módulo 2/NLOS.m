% Channel gain NLOS
clear all
close all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Room Specs
Xmax = 6;               % Max width
Ymax = 6;               % Max length
h = 3;                  % Ceiling height
rho = 0.7;              % wall reflection coefficient«
dx = 0.06;              % spatial resolution for ground dy=dx
dr = 0.06;              % spatial resolution for reflections dy_r=dy
A_r = dr^2;             % area elements for reflections
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TX and RX Specs
%m = - log(2)/log(cosd(60));
m = 1;
FOV = pi/2;
R = 0.7;               % RX responsivity
Apd = 100e-6;           % area of the photo-detector
Pt = 1;                % transmitted power
% TX positions
[XT, YT] = meshgrid([-2 2],[-2 2]);
NT = sqrt(prod(size(XT)));
nt = [0 0 -1];

% RX positions
[X, Y] = meshgrid(-Xmax/2+dx/2:dx:Xmax/2-dx/2,-Ymax/2+dx/2:dx:Ymax/2-dx/2);
NR = sqrt(prod(size(X)));
nr = [0 0 1];

% NLOS contribution
Pr_NLOS_TX = cell(NT,NT);


CN_LOS = rho * R * Pt * (m+1)/(2*pi*pi) * A_r*Apd;

nw = [1 0 0];

Xw = -Xmax/2+dx/2:dx:Xmax/2-dx/2;
Yw = -Ymax/2+dx/2:dx:Ymax/2-dx/2;
Zw = dr/2:dr:h-dr/2;


Nwxyz = numel(Xw);
Nwz = numel(Zw);


for t1 = 1:NT
    for t2 = 1:NT
        Pr = zeros(size(X));
        Tx_pos = [XT(t1,t2) YT(t1,t2) h];
        for r1 = 1:NR
            for r2 = 1:NR
                Rx_Pos = [X(r1,r2) Y(r1,r2) 0];
                for w1 = 1:Nwxyz
                    for w2=1:Nwz

                        Wx = [1 0 0];
                        W_pos = [-Xmax/2 Yw(w1) Zw(w2)];

                        d1 = [W_pos(1)-Tx_pos(1), W_pos(2)-Tx_pos(2), 
                             W_pos(3)-Tx_pos(3)];
                        d2 = [Rx_Pos(1) - W_pos(1), Rx_Pos(2) - W_pos(2),
                             Rx_Pos(3) - W_pos(3)];
    
                        d1n = sqrt(d1*d1.');
                        d2n = sqrt(d2*d2.');
    
                        cos_phi = (nt*d1.')/d1n;
                        cos_alpha = -(nw*d1.')/d1n;
                        cos_beta = (nw*d2.')/d2n;
                        cos_psi = -(nr*d2.')/d2n;
        
                        if acos(cos_psi) < FOV
                            Pr(r1,r2) = Pr(r1,r2) + CN_LOS*((cos_phi^m)*cos_alpha)
                                        /((d1n)^2) * (cos_beta*cos_psi)/(d2n^2);
                        end
                    end
                end
            end
        end
        Pr_NLOS_TX{t1,t2} = Pr;
    end
end

figure(1)
Pr_NLOS = PlotChannel(X,Y,NT,Pr_NLOS_TX,'P_r_N_L_O_S_1 (W)');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Auxiliary functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function final = PlotChannel(X,Y,NT,C,str)
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

    final = zeros(length(X),length(Y));
    for i=1:4
        final = final + rot90(M,i);
    end

    surfc(X,Y,final)
    shading interp
    xlabel('x (m)')
    ylabel('y (m)')
    zlabel(str)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%