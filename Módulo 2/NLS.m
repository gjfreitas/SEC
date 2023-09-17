%non linear least squares position estimation
clear all
close all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load channel_6cm.mat
% Room Specs
Xmax = 6;               % Max width
Ymax = 6;               % Max length
h = 3;                  % Ceiling height
% TX and RX Specs
hpa = 60;               % TX half power angle
m = - log(2)/log(cosd(hpa));
FOV = 90;
R = 0.7;                % RX responsivity
Apd = 100e-6;           % area of the photo-detector
Pt = 1;                 % transmitted power
SNR = 240;               % Signal to noise ratio

N = 4; %graus do polinomio de interpolacao
NT = sqrt(NTx);

Navg = 10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Pr_TX = cell(sqrt(NTx),sqrt(NTx));
d_TX = cell(sqrt(NTx),sqrt(NTx));
r = cell(sqrt(NTx),sqrt(NTx));

for p = 1:Navg % averaging
    NoiseF=1+(10^(-SNR/20))*randn(size(Pr_LOS_TX{1,1}));
    for t1 = 1:sqrt(NTx)
        for t2 = 1:sqrt(NTx)
            Pr_TX{t1,t2} = (Pr_LOS_TX{t1,t2} + Pr_NLOS1_TX{t1,t2}).* NoiseF;
            coeff = DistanceEstimator(Pr_TX,XT,YT,NT,X,Y,N,SNR);
            d_TX{t1,t2} = polyval(coeff,Pr_TX{t1,t2});
            r{t1,t2} = d_TX{t1,t2}.^2-h.^2;
        end
    end
end

X_ = zeros(size(X));
Y_ = zeros(size(Y));
for r1 = 1:sqrt(numel(X))
    for r2 = 1:sqrt(numel(Y))
        A = [XT(1,2) - XT(1,1) YT(1,2) - YT(1,1); XT(2,1) - XT(1,1) 
            YT(2,1) - YT(1,1); XT(2,2) - XT(1,1) YT(2,2) - YT(1,1)];
        B = [r{1,1}(r1,r2)^2 - r{1,2}(r1,r2)^2 + XT(1,2)^2 + YT(1,2)^2
            - XT(1,1)^2 - YT(1,1)^2];
        B = [B; r{1,1}(r1,r2)^2 - r{2,1}(r1,r2)^2 + XT(2,1)^2 + YT(2,1)^2
            - XT(1,1)^2 - YT(1,1)^2];
        B = 0.5*[B; r{1,1}(r1,r2)^2 - r{2,2}(r1,r2)^2 + XT(2,2)^2 + YT(2,2)^2
            - XT(1,1)^2 - YT(1,1)^2];
        
        x0 = pinv(A)*B;
        X0_(r1,r2) = x0(1);
        Y0_(r1,r2) = x0(2);
        options = optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective');
        fun = @(x)Jcost3D(x,XT,YT,NT,h,d_TX,r1,r2);
        x = lsqnonlin(fun,[x0;0],[],[],options);
        X_(r1,r2) = x(1);
        Y_(r1,r2) = x(2);
    end
end
% Error surface
err_Surf = sqrt((X-X_).^2 + (Y-Y_).^2);
figure(1)
surf(X,Y,err_Surf)
xlabel('x')
ylabel('y')
zlabel('Err')
% Cumulative distribution function of the error
figure(2)
[cnt,x] = hist(err_Surf(:),100);
cdf = cumsum(cnt/numel(err_Surf));
plot(x,cdf)
xlabel('err')
ylabel('CDF')
grid on 

function p = DistanceEstimator(Pr_TX,XT,YT,NT,X,Y,N,SNR)
    Pr = Pr_TX{1,1}.*(1+(10^(-SNR/20))*randn(size(X)));
    d = sqrt((X-XT(1,1)).^2+(Y-YT(1,1)).^2);
    p = polyfit(Pr,d,N);
end

function c = Jcost3D(x,XT,YT,NT,h,d_TX,r1,r2)
    xt = reshape(XT,numel(XT),1);
    yt = reshape(YT,numel(YT),1);
    dTX = zeros(size(XT));
    for t1=1:NT
        for t2=1:NT
            dTX(t1,t2)=d_TX{t1,t2}(r1,r2);
        end
    end
    dtx = reshape(dTX,numel(dTX),1);
    c = sqrt((x(1)-xt).^2+(x(2)-yt).^2+(x(3)-h).^2)-dtx;
end


function qerr = QauntileCDF(err,Q)
    N = numel(err);
    [cnt,x] = hist(err(:),100);
    cdf = cumsum(cnt/N);
    flag = 1;
    eps = 0.005;

    while(flag==1)
        if isempty(x(abs(cdf-Q) < eps))
            eps = eps + 0.005;
        else
            qerr = x(abs(cdf-Q)<eps);
            flag = 0;
        end
    end
end