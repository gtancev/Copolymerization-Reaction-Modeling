%% Polymerization Reaction Model
% George Tancev

clear all; close all; clc;

%% data

V = 15;                     % L
cA_0 = 3.5;                 % mol/L
pA = 0.90;                  % kg/L
pB = 0.94;                  % kg/L
wI_0 = 0.01;    
MA = 0.104;                 % kg/mol
MB = 0.1;                   % kg/mol
MI = 0.164;                 % kg/mol

mA_0 = cA_0*V*MA;           % kg
VA_0 = mA_0/pA;             % L
mB_0 = (V-VA_0)*pB;         % kg
cB_0 = mB_0/(MB*V);         % mol/L

X_A_0 = cA_0/(cA_0+cB_0);

mI_0 = wI_0*(mA_0+mB_0);    % kg
cI_0 = mI_0/(MI*V);         % mol/L

f = 0.5;
C_eta = 25;
C_RD = 180;

kd = 6.77e-6;               % s^(-1)

kp_AA_0 = 4.1e2;            % L/(mol s)
kp_AD_0 = 2e11;             % L/(mol s)
kt_AA_0 = 2.4e7;            % L/(mol s)
kt_AD_0 = 5e8;              % L/(mol s)

kp_BB_0 = 9.3e2;            % L/(mol s)
kp_BD_0 = 2e11;             % L/(mol s)
kt_BB_0 = 9.2e6;            % L/(mol s)
kt_BD_0 = 5e8;              % L/(mol s)

rA = 0.52;
rB = 0.46;

kp_AB_0 = kp_AA_0/rA;       % L/(mol s)
kp_BA_0 = kp_BB_0/rB;       % L/(mol s)

%% Solving the system of ODEs

tspan = 1:1:7200;           % s
[t,c] = ode15s(@(t,c)batch(t,c),tspan,[cI_0 cA_0 0 cB_0 0]);
t = t/3600;

%% a)

X_A = linspace(0,1,100);
F_A = ((rA-1).*X_A.^2+X_A)./((rA-2).*X_A.^2+2.*X_A+rB.*(1-X_A).^2);
DIAG = X_A;

x = 1-(c(:,2)+c(:,4))/(cA_0+cB_0);

F_A_c = cumulative( t,c );

figure(1);
subplot(3,2,1);
plot(X_A,F_A,X_A,DIAG,'--');
title('Mayo-Lewis diagram');
xlabel({'X_A'});
ylabel({'F_A'});

figure(1);
subplot(3,2,2);
plot(x(1:end-1),F_A_c(1:end,1),x(3:end),F_A_c(2:end,2));
axis([0 1 0 0.6])
title('cumulative composition distribution');
xlabel({'conversion','x'});
ylabel({'F_A, F_A^c'});
legend('F_A','F_A^c','Location','best');

%% b)

figure(1);
subplot(3,2,3);
plot(t,c(:,2),t,c(:,4));
axis([0 2 0 6]);
title('concentration profile');
xlabel({'time','h'});
ylabel({'concentration','mol / L'});
legend('styrene','methyl methacrylate','Location','best');

figure(1);
subplot(3,2,4);
plot(t,x);
title('conversion profile');
xlabel({'time','h'});
ylabel({'conversion','x'});

%% c)

figure(1);
subplot(3,2,5);
plot(x(3:end),F_A_c(2:end,3),x(3:end),F_A_c(2:end,4),x(3:end),F_A_c(2:end,6),x(3:end),F_A_c(2:end,7));
title({'cumulative and instantanenous', 'number/ weight average'});
xlabel({'conversion','x'});
ylabel({'n_N, n_W'});
legend('n^c_N','n^c_W','n_N','n_W','Location','best');

figure(1);
subplot(3,2,6);
plot(x(3:end),F_A_c(2:end,5),x(3:end),F_A_c(2:end,8));
title({'cumulative and instantaneous', 'polydispersity'});
xlabel({'conversion','x'});
ylabel({'\sigma'});
legend('\sigma^c','\sigma','Location','best');

%%
