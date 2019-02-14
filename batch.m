function [ F ] = batch( t,c )

cI = c(1);
cA = c(2);
cA_R = c(3);
cB = c(4);
cB_R = c(5);

%%

V = 15;                     % L
cA_0 = 3.5;                 % mol/L
pA = 0.90;                  % kg/L
pB = 0.94;                  % kg/L
MA = 0.104;                 % kg/mol
MB = 0.1;                   % kg/mol

mA_0 = cA_0*V*MA;           % kg
VA_0 = mA_0/pA;             % L
mB_0 = (V-VA_0)*pB;         % kg
cB_0 = mB_0/(MB*V);         % mol/L

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

%%

wp = (MA*V*(cA_0-c(2))+MB*V*(cB_0-c(4)))/(mA_0+mB_0);

kp_AA = 1/(1/kp_AA_0 + exp(C_eta*wp)/kp_AD_0);
kt_AA = 1/(1/kt_AA_0 + exp(C_eta*wp)/kt_AD_0) + C_RD*kp_AA*(1-wp);

kp_BB = 1/(1/kp_BB_0 + exp(C_eta*wp)/kp_BD_0);
kt_BB = 1/(1/kt_BB_0 + exp(C_eta*wp)/kt_BD_0) + C_RD*kp_BB*(1-wp);

kp_AB = kp_AA/rA;
kt_AB = sqrt(kt_AA*kt_BB);

kp_BA = kp_BB/rB;
kt_BA = kt_AB;

%%

dcI = - kd*cI;
dcA = - kp_AA*cA*cA_R - kp_BA*cA*cB_R;
dcA_R = + 2*f*kd*cI*(cA/(cA+cB)) - kt_AA*cA_R^2 - kt_AB*cA_R*cB_R - kp_AB*cA_R*cB + kp_BA*cB_R*cA;
dcB = - kp_BB*cB*cB_R - kp_AB*cB*cA_R;
dcB_R = + 2*f*kd*cI*(cB/(cA+cB)) - kt_BB*cB_R^2 - kt_BA*cA_R*cB_R - kp_BA*cB_R*cA + kp_AB*cA_R*cB;

%%

F = [dcI; dcA; dcA_R; dcB; dcB_R];

end

