function [ D ] = cumulative( t,c )

%% data

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
    
    C_t = 1000;
    
%%

    h = length(t);

    for t1=1:1:h;

    c_A = c(t1,2);
    c_A_dot = c(t1,3);
    c_B = c(t1,4);
    c_B_dot = c(t1,5);
    c_R_dot = c_A_dot+c_B_dot;
    
    wp = (MA*V*(cA_0-c_A)+MB*V*(cB_0-c_B))/(mA_0+mB_0);
    
    kp_AA = 1/(1/kp_AA_0 + exp(C_eta*wp)/kp_AD_0);
    kt_AA = 1/(1/kt_AA_0 + exp(C_eta*wp)/kt_AD_0) + C_RD*kp_AA*(1-wp);

    kp_BB = 1/(1/kp_BB_0 + exp(C_eta*wp)/kp_BD_0);
    kt_BB = 1/(1/kt_BB_0 + exp(C_eta*wp)/kt_BD_0) + C_RD*kp_BB*(1-wp);

    kp_AB = kp_AA/rA;
    kt_AB = sqrt(kt_AA*kt_BB);

    kp_BA = kp_BB/rB;
    kt_BA = kt_AB;
    
    p_A = c_A_dot/c_R_dot;
    p_B = c_B_dot/c_R_dot;
    
    kpM = (kp_AA*p_A+kp_BA*p_B)*c_A+(kp_BB*p_B+kp_AB*p_A)*c_B;
    kt = kt_AA*p_A^2+2*kt_AB*p_A*p_B+kt_BB*p_B^2;
    
    ktc = kt/(1+C_t);
    ktd = kt-ktc;
    
    beta = ktc*c_R_dot/(kpM);
    gamma = ktd*c_R_dot/(kpM);
    alph = beta + gamma;
    
    R_p_A = (kp_AA*p_A+kp_BA*p_B)*c_A*c_R_dot;
    R_p_B = (kp_BB*p_B+kp_AB*p_A)*c_B*c_R_dot;
    R_p = R_p_A+R_p_B;
    
    dPdt(t1) = R_p*(gamma+beta/2);
    
    n_N(t1) = 1/(gamma+beta/2);
    n_W(t1) = 2*(gamma+1.5*beta)/alph^2;
    S1(t1) = n_W(t1)/n_N(t1);
    
    mu1(t1) = n_N(t1)*dPdt(t1);
    mu2(t1) = n_W(t1)*mu1(t1);
    
    F_A_1(t1) = (((rA*c_A+c_B)*c_A)/((rA*c_A+c_B)*c_A+(rB*c_B+c_A)*c_B));
    F_A(t1) = (((rA*c_A+c_B)*c_A)/((rA*c_A+c_B)*c_A+(rB*c_B+c_A)*c_B))*dPdt(t1);
    
    end
    
    tspan = 0:1:length(t);

    for t1=2:1:length(t)-1;
        
        P_t(t1) = trapz(tspan(2:t1+1),dPdt(2:t1+1));
        F_A_c(t1) = trapz(tspan(2:t1+1),F_A(2:t1+1))/P_t(t1);
        mu1c(t1) = trapz(tspan(2:t1+1),mu1(2:t1+1))/P_t(t1);
        mu2c(t1) = trapz(tspan(2:t1+1),mu2(2:t1+1))/P_t(t1);
        n_N_c(t1) = mu1c(t1);
        n_W_c(t1) = mu2c(t1)/mu1c(t1);
        S2(t1) = mu2c(t1)/(mu1c(t1)^2);

    end
    
    D = [F_A_1(1:end-1)' F_A_c' n_N_c' n_W_c' S2' n_N(1:end-1)' n_W(1:end-1)' S1(1:end-1)'];


end