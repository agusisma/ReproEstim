%% Code by Agus Hasan
% Estimation of the time-varying reproduction number using active cases
% based on SIR model

%%
clear;
clc;

%% Inputs
% primary inputs (you can change this input according to your data)
load DATADE.txt;                                      % data structure: date | month | active case | removed case (recovered+deceased)
DATA = DATADE;
%N        = 5792202;                                % number of population Denmark
N        = 83783942;                                % number of population Germany

load CORIDE.txt;

interval = CORIDE;

Tinf     = 12;                                      % infectious time
Std_Tinf = 1;                                       % standard deviation of infectious time
CFR      = 0.04;                                    % case fatality rate
% secondary inputs (you don't have to change this input)
tf       = length(DATA);                            % simulation time
td       = datetime(2020,DATA(1,2),DATA(1,1)-1) + caldays(1:tf); % time horizon
dt       = 0.01;                                    % time steps
t        = dt:dt:tf;
sigma    = 1.96;                                    % CI
C        = [1 0 0 0;                                % data matrix
            0 1 0 0;
            0 0 1 0];
QF       = diag([10 10 10 0.2]);                    % tuning gain 1
RF       = diag([100 10 1]);                        % tuning gain 2
% calculate susceptible population
for i=1:tf
    DATA(i,5) = N-sum(DATA(i,3:4));
end


%L = [0.3025 -0.084 0; %DK
%    -0.0084 0.64   0;
%    0       0      0.916;
%    -0.023  0.064  0];

L = [0.405  -0.6   0; %DE
    -0.06   0.94   0;
    0       0      0.92;
    -0.0004 0.0032 0];


%% For plotting
% Adding a low pass-filter to handle short-term data fluctuations 

windowSize = 700;
b          = (1/windowSize)*ones(1,windowSize);
a          = 1;
% Plotting Rt below 1
curve11 = 0*ones(1,tf);
curve22 = 1*ones(1,tf);
x2      = [td, fliplr(td)];
fs      = 48;

%% Simulation
for j = 1:3
% Infection time
Ti     = Tinf-Std_Tinf*sigma+(j-1)*Std_Tinf*sigma;
gamma  = 1/Ti;

%% Initialization
xhat     = [N-1; 1; 0; 1];   % initial condition
xbar     = [N-1; 1; 0; 1];   % initial condition
Pplus    = 1000*eye(4);            % since we know excatly the initial conditions
xhatEff  = 0;
xbarEff  = 0;
% for plotting
xArray       = [];
xhatArray    = [];
xbarArray    = [];
xhatEffArray = [];
xbarEffArray = [];
% extended Kalman filter
for i=1:((tf-1)/dt)
    xhatArray    = [xhatArray xhat]; 
    xhatEffArray = [xhatEffArray xhatEff];          
    xbarArray    = [xbarArray xbar]; 
    xbarEffArray = [xbarEffArray xbarEff];      
    % assimilating the reported data
    y = [interp1(0:1:tf-1,DATA(:,5),t);
        interp1(0:1:tf-1,DATA(:,3),t);
        interp1(0:1:tf-1,DATA(:,4),t)];
    % prediction
    xhat(1) = xhat(1)-gamma*xhat(4)*xhat(1)*xhat(2)*dt/N;
    xhat(2) = xhat(2)+gamma*xhat(4)*xhat(1)*xhat(2)*dt/N-gamma*xhat(2)*dt;
    xhat(3) = xhat(3)+gamma*xhat(2)*dt;
    xhat(4) = xhat(4);
    % calculating the Jacobian matrix
    FX    = [1-gamma*xhat(4)*xhat(2)*dt/N -gamma*xhat(4)*xhat(1)*dt/N 0 -gamma*xhat(1)*xhat(2)*dt/N;
             gamma*xhat(4)*xhat(2)*dt/N 1+gamma*xhat(4)*xhat(1)*dt/N-gamma*dt 0 gamma*xhat(1)*xhat(2)*dt/N;
             0 gamma*dt 1 0;
             0 0 0 1];
    Pmin  = FX*Pplus*FX'+QF;
    % update 
    KF    = Pmin*C'*inv(C*Pmin*C'+RF);
    xhat  = xhat + KF*(y(:,i)-C*xhat);
    Pplus = (eye(4)-KF*C)*Pmin;
    xhat(4) = max(0,xhat(4));
    xhatEff = (xhat(1)/N)*xhat(4);

    % NLO
    
    xbar(1) = xbar(1)-gamma*xbar(4)*xbar(1)*xbar(2)*dt/N;
    xbar(2) = xbar(2)+gamma*xbar(4)*xbar(1)*xbar(2)*dt/N-gamma*xbar(2)*dt;
    xbar(3) = xbar(3)+gamma*xbar(2)*dt;
    xbar(4) = xbar(4);
    
    xbar  = xbar + L*(y(:,i)-C*xbar);
    xbar(4) = max(0,xbar(4));
    xbarEff = (xbar(1)/N)*xbar(4);
    
end

%% Plotting

xhatArray(4,:) = filter(b,a,xhatEffArray);

xhatSArray  = [];
xhatS       = xhatArray(1,tf);
xhatIArray  = [];
xhatI       = xhatArray(2,tf);
xhatRArray  = [];
xhatR       = xhatArray(3,tf);
xhatRtArray = [];
xhatRt      = xhatArray(4,tf);
for i=1:tf-1
    xhatSArray  = [xhatSArray xhatS];
    xhatS       = xhatArray(1,(1/dt)*i);
    xhatIArray  = [xhatIArray xhatI];
    xhatI       = xhatArray(2,(1/dt)*i);
    xhatRArray  = [xhatRArray xhatR];
    xhatR       = xhatArray(3,(1/dt)*i);
    xhatRtArray = [xhatRtArray xhatRt];
    xhatRt      = xhatArray(4,(1/dt)*i);
end

xhatSArray  = [xhatSArray xhatS];
xhatIArray  = [xhatIArray xhatI];
xhatRArray  = [xhatRArray xhatR];
xhatRtArray = [xhatRtArray xhatRt];

M(j,:) = xhatRtArray;

xbarArray(4,:) = filter(b,a,xbarEffArray);

xbarSArray  = [];
xbarS       = xbarArray(1,tf);
xbarIArray  = [];
xbarI       = xbarArray(2,tf);
xbarRArray  = [];
xbarR       = xbarArray(3,tf);
xbarRtArray = [];
xbarRt      = xbarArray(4,tf);
for i=1:tf-1
    xbarSArray  = [xbarSArray xbarS];
    xbarS       = xbarArray(1,(1/dt)*i);
    xbarIArray  = [xbarIArray xbarI];
    xbarI       = xbarArray(2,(1/dt)*i);
    xbarRArray  = [xbarRArray xbarR];
    xbarR       = xbarArray(3,(1/dt)*i);
    xbarRtArray = [xbarRtArray xbarRt];
    xbarRt      = xbarArray(4,(1/dt)*i);
end

xbarSArray  = [xbarSArray xbarS];
xbarIArray  = [xbarIArray xbarI];
xbarRArray  = [xbarRArray xbarR];
xbarRtArray = [xbarRtArray xbarRt];

V(j,:) = xbarRtArray;

end

for l = 1:tf
    curve2(l)      = max(M(:,l));
    xhatRtArray(l) = mean(M(:,l));
    curve1(l)      = min(M(:,l));
end

for l = 1:tf
    curve2nlo(l)      = max(V(:,l));
    xbarRtArray(l)    = mean(V(:,l));
    curve1nlo(l)      = min(V(:,l));
end

%% Plotting
% 
% figure(1)
% subplot(2,1,1)
% plot(td,DATA(:,3),'ob','LineWidth',6)
% hold on
% plot(td,xhatIArray,'-r','LineWidth',4);
% xlim([min(td)+7 max(td)])
% ylabel('Active Case')
% set(gca,'color','none','FontSize',36)
% legend('Reported','Estimated','Location','northwest')
% grid on;
% grid minor;
% subplot(2,1,2)
% plot(td,DATA(:,4),'ob','LineWidth',6);
% hold on
% plot(td,xhatRArray,'-r','LineWidth',4)
% xlim([min(td)+7 max(td)])
% ylabel('Removed Case')
% set(gca,'color','none','FontSize',36)
% legend('Reported','Estimated','Location','northwest')
% grid on;
% grid minor;

figure(2)
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, 'r');
alpha(0.5)
hold on;
plot(td,xhatRtArray,'-r','LineWidth',4)
hold on
inBetween = [curve1nlo, fliplr(curve2nlo)];
fill(x2, inBetween, 'g');
alpha(0.5)
hold on;
plot(td,xbarRtArray,'g','LineWidth',4)
hold on;
plot(td,ones(1,tf),':b','LineWidth',4)
set(gca,'color','none','FontSize',36)
%ylim([0.5 2])
xlim([min(td)+15 max(td)-15])
ylabel('Rt','color','k')
grid on
grid minor
legend('EKF','mean EKF','NLO','mean NLO')