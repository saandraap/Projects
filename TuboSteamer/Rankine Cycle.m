%% RANKINE DUAL

% UNIDADES acordes con XSteam function:

% m     Flujo masico (kg/s)
% T     Temperature (°C)
% p	    Pressure    (bar)
% h	    Enthalpy    (kJ/kg)
% v	    Specific volume  (m3/kg)
% rho	Density  (kg/ m3) 
% s	    Specific entropy  (kJ/ kg K)
% u	    Specific internal energy  (kJ/kg)
% Cp	Specific isobaric heat capacity  (kJ/ kg K)
% Cv	Specific isochoric heat capacity  (kJ/ kg K)
% my	Viscosity   (Pa s)
% tc	Thermal Conductivity ( W/m °C)
% x	    Vapour fraction
% vx	Vapour Volume Fraction

% FACTOR DE CONVERSION 

bar_Pa = 100000;
Pa_bar = 0.00001;
w_hp = 0.00134102;

%% gases de escape del motor

% designado como: e

m_e = 0.054;
TeK = 720;
Te = TeK-273;
cp_mathcad = 1.241;
cp_e= airProp(Te, 'cp');

%% fluido refrigerante 

% agua, designado como: ref 

m_ref = 0.4;
Tref = 100;
p_ref = 1;

%% CICLO ALTA TEMPERATURA

% condiciones iniciales, liquido saturado.
% fluido de trabajo: agua, designado como: a

m_a = 0.03;
p1 = 1;
n = 0;

%% CICLO BAJA TEMPERATURA

% condiciones iniciales, liquido saturado.
% fluido de trabajo: agua, designado como: w

m_w = 0.03; 
T5 = 25;


%% propiedades estado 1

T1 = XSteam ('Tsat_p', p1);
h1 = XSteam ('hL_p', p1);
s1 = XSteam ('sL_p', p1);
cp1 = XSteam ('cpL_p', p1);

%% propiedades estado 2 - Incremento de presion isoentropico hasta 7 bar.

p2 = 7;
s2_id = s1;
h2_id = XSteam('h_ps', p2, s2_id);
Wb1_id = m_a*(h2_id-h1);       % kW
n_bomb = 0.45; % ******************* COMPROBAR **********

Wb1_real = Wb1_id/n_bomb;
h2 = h1 + (Wb1_real/m_a);
T2 = XSteam ('T_ph', p2, h2);
s2 = XSteam ('s_ph', p2,h2);
Wb1_real_hp = Wb1_real*1000*w_hp;

%% EVAPORADOR ALTA TEMPERATURA - INTERCAMBIADOR DE PLACAS
%  Intercambiador numero 1

Np = 25;
Lp = 0.25;
Wp = 0.15;
tp = 0.0004;
bp = 0.002;
V = 1.10;
At = Np*Lp*Wp;
De = 2*bp/V;    % Igual para agua y para gases de escape.

% AGUA

Te_a1 = (((Te+T2)/2)+T2)/2;  % Temperatura equivalente para propiedades del agua
cp_a1 = XSteam ('Cp_pt', p2, Te_a1);
my_a1 = XSteam ('my_pt', p2, Te_a1);
k_a1 = XSteam ('tc_pt', p2, Te_a1);
Re_a1 = 2*m_a*De / (bp*Wp*Np*my_a1);
Pr_a1 = cp_a1*my_a1/(k_a1/1000);

Nu_a1 = 0.2946*Re_a1^0.7*Pr_a1^(1/3);
h_a1 = Nu_a1*k_a1/De;    % (W/m2 K)

% GAS

Te_g1 = (((Te+T2)/2)+Te)/2;   % Temperatura equivalente para propiedades del gas
cp_g1 = airProp (Te_g1, 'cp');
my_g1 = airProp2 (Te_g1, 'my');
k_g1 =  airProp2 (Te_g1, 'k');
Re_g1 = 2*m_e*De / (bp*Wp*Np*my_g1);
Pr_g1 = cp_g1*my_g1/(k_g1/1000);

Nu_g1 = 0.2946*Re_g1^0.7*Pr_g1^(1/3);
h_g1 = Nu_g1*k_g1/De;    % (W/m2 K)


% Calculo intercambiador

UA1 = 1/(1/(h_a1*At)+1/(h_g1*At));   %(W/K)

C_a1 = m_a*cp_a1*10^3;
C_g1 = m_e*cp_g1*10^3;
C_min1 = min(C_a1,C_g1);
C_max1 = max(C_a1,C_g1);
C_r1 = C_min1/C_max1;
C_r_cond = 0;
NTU1 = UA1/C_min1;
Efic1 = (1-exp(-NTU1*(1-C_r1)))/(1-C_r1*exp(-NTU1*(1-C_r1)));
Qint1 = Efic1*C_min1*(Te-T2);   % (W)    %***************
T3 = T2 + Qint1*10^-3/(m_a*cp_a1);

% Perdida de presion

C_fric_a1 = 0.37*Re_a1^-0.172;
G_a1 = 2*m_a / (bp*Wp*Np);
dens_a1 = XSteam ('rho_pT', p2, Te_a1);
Perd_pres1 = 4*C_fric_a1*G_a1^2*Lp/(2*dens_a1*De);
Perd_pres_bar1 = Perd_pres1 * Pa_bar;
p3 = p2-Perd_pres_bar1;

%% PROPIEDADES ESTADO 3

h3 = XSteam ('h_PT', p3, T3);
s3 = XSteam ('s_PT', p3, T3);
x3 = XSteam ('x_ph', p3, h3);

%% propiedades estado 4 - Expansion  ( Pi_tur = 7)

p4 = p1;
s4_id = s3;
h4_id = XSteam ('h_ps', p4, s4_id);
Wtur1_id = m_a*(h3-h4_id);

n_tur = 0.45; %***** variar con T********
Wtur1_real = Wtur1_id*n_tur; 
Wtur1_real_hp = Wtur1_real*1000*w_hp;

h4 = h3-(Wtur1_real/m_a);
s4 = XSteam ('s_ph', p4, h4);
T4 = XSteam ('T_ps', p4, s4);
x4 = XSteam ('x_ph', p4, h4);

incr_p = p3-p4;

%% Propiedades estado 5 

p5 = XSteam ('Psat_T',T5);
h5 = XSteam ('hL_p', p5);
s5 = XSteam ('sL_p', p5);

%% Propiedades estado 6 - Incremento de presion isoentropico 

p6 = XSteam ('Psat_T', 80);
s6_id = s5;
h6_id = XSteam('h_ps', p6, s6_id);
Wb2_id = m_w*(h6_id-h5);       % kW
n_bomb2 = 0.45; % ******************* COMPROBAR **********

Wb2_real = Wb2_id/n_bomb2;
h6 = h5 + (Wb2_real/m_w);
T6 = XSteam ('T_ph', p6, h6);
s6 = XSteam ('s_ph', p6,h6);
Wb2_real_hp = Wb2_real*1000*w_hp;

%% REFRIGERACION- INTERCAMBIADOR DE PLACAS (6-7)
%  Intercambiador numero 2

Np = 25;
Lp = 0.25;
Wp = 0.15;
tp = 0.0004;
bp = 0.002;
V = 1.10;
At = Np*Lp*Wp;
De = 2*bp/V;    % Igual para ambos lados

% AGUA CICLO

Te_w2 = (((Tref+T6)/2)+T6)/2;  % Temperatura equivalente para propiedades del agua
cp_w2 = XSteam('Cp_pt', p6, Te_w2);
my_w2 = XSteam('my_pt', p6, Te_w2);
k_w2 = XSteam('tc_pt', p6, Te_w2);
Re_w2 = 2*m_w*De / (bp*Wp*Np*my_w2);
Pr_w2 = cp_w2*my_w2/(k_w2/1000);

Nu_w2 = 0.2946*Re_w2^0.7*Pr_w2^(1/3);
h_w2 = Nu_w2*k_w2/De;    % (W/m2 K)

% REFRIGERANTE  (AGUA)

Te_ref2 = (((Tref+T2)/2)+Tref)/2; % Temperatura equivalente para propiedades del agua de refrigeracion
cp_ref2 = XSteam('Cp_pt', p_ref, Te_ref2);
my_ref2 = XSteam('my_pt', p_ref, Te_ref2);
k_ref2 = XSteam('tc_pt', p_ref, Te_ref2);

Re_ref2 = 2*m_ref*De / (bp*Wp*Np*my_ref2);
Pr_ref2 = cp_ref2*my_ref2/(k_ref2/1000);

Nu_ref2 = 0.2946*Re_ref2^0.7*Pr_ref2^(1/3);
h_ref2 = Nu_ref2*k_ref2/De;    % (W/m2 K)


% Calculo intercambiador

UA2 = 1/(1/(h_ref2*At)+1/(h_w2*At));   %(W/K)

C_w2 = m_w*cp_w2*10^3;
C_ref2 = m_ref*cp_ref2*10^3;
C_min2 = min(C_w2,C_ref2);
C_max2 = max(C_w2,C_ref2);
C_r2 = C_min2/C_max2;
NTU2 = UA2/C_min2;
Efic2 = (1-exp(-NTU2*(1-C_r2)))/(1-C_r2*exp(-NTU2*(1-C_r2)));
Qint2 = Efic2*C_min2*(Tref-T6);   % (W)    %***************
T7 = T6 + Qint2*10^-3/(m_w*cp_w2);

% Perdida de presion

C_fric_w2 = 0.37*Re_w2^-0.172;
G_w2 = 2*m_w / (bp*Wp*Np);
dens_w2 = XSteam ('rho_pT', p6, Te_w2);
Perd_pres2 = 4*C_fric_w2*G_w2^2*Lp/(2*dens_w2*De);
Perd_pres_bar2 = Perd_pres2 * Pa_bar;
p7 = p6-Perd_pres_bar2;

%% Propiedades estado 7

h7 = XSteam ('h_PT', p7, T7);
s7 = XSteam ('s_PT', p7, T7);
x7 = XSteam ('x_ph', p7, T7);

%% Evaporador ciclo baja temperatura (7-8) y (4-4') con calor residual ciclo de alta.
%  Intercambiador numero 3

Np = 25;
Lp = 0.25;
Wp = 0.15;
tp = 0.0004;
bp = 0.002;
V = 1.10;
At = Np*Lp*Wp;
De = 2*bp/V;    % Igual para ambos lados

% AGUA CICLO BAJA TEMP.

Te_w3 = (((T4+T7)/2)+T7)/2;  % Temperatura equivalente para propiedades del agua
cp_w3 = XSteam('Cp_pt', p6, Te_w3);
my_w3 = XSteam('my_pt', p6, Te_w3);
k_w3 = XSteam('tc_pt', p6, Te_w3);
Re_w3 = 2*m_w*De / (bp*Wp*Np*my_w3);
Pr_w3 = cp_w3*my_w3/(k_w3/1000);

Nu_w3 = 0.2946*Re_w3^0.7*Pr_w3^(1/3);
h_w3 = Nu_w3*k_w3/De;    % (W/m2 K)


% AGUA CICLO ALTA TEMP.

Te_a3 = (((T4+T7)/2)+T4)/2;  % Temperatura equivalente para propiedades del agua
cp_a3 = XSteam ('Cp_pt', p2, Te_a3);
my_a3 = XSteam ('my_pt', p2, Te_a3);
k_a3 = XSteam ('tc_pt', p2, Te_a3);
Re_a3 = 2*m_a*De / (bp*Wp*Np*my_a3);
Pr_a3 = cp_a3*my_a3/(k_a3/1000);

Nu_a3 = 0.2946*Re_a3^0.7*Pr_a3^(1/3);
h_a3 = Nu_a3*k_a3/De;    % (W/m2 K)


% Calculo intercambiador

UA3 = 1/(1/(h_a3*At)+1/(h_w3*At));   %(W/K)

C_a3 = m_a*cp_a3*10^3;
C_w3 = m_w*cp_w3*10^3;
C_min3 = min(C_a3,C_w3);
C_max3 = max(C_a3,C_w3);
C_r3 = C_min3/C_max3;
C_r_ev = 0;
NTU3 = UA3/C_min3;
Efic3 = (1-exp(-NTU1*(1-C_r_ev)))/(1-C_r_ev*exp(-NTU1*(1-C_r_ev)));
Qint3 = Efic3*C_min3*(T4-T7);   % (W)    %***************
T8 = T7 + Qint3*10^-3/(m_w*cp_w3);
T4b = T4 - Qint3*10^-3/(m_a*cp_a3);


% Perdida de presion

C_fric_w3 = 0.37*Re_w3^-0.172;
G_w3 = 2*m_w / (bp*Wp*Np);
dens_w3 = XSteam ('rho_pT', p7, Te_w3);
Perd_pres3 = 4*C_fric_w3*G_w3^2*Lp/(2*dens_w3*De);
Perd_pres_bar3 = Perd_pres3 * Pa_bar;
p8 = p7-Perd_pres_bar3;


%% PROPIEDADES ESTADO 4'

p4b = p4;
h4b = XSteam ('h_PT', p4b, T4b);
s4b = XSteam ('s_PT', p4b, T4b);
x4b = XSteam ('x_ph', p4b, T4b);

%% PROPIEDADES ESTADO 8

h8 = XSteam ('h_PT', p8, T8);
s8 = XSteam ('s_PT', p8, T8);
x8 = XSteam ('x_ph', p8, h8);

%% propiedades estado 9 - Expansion  

p9 = p5;
s9_id = s8;
h9_id = XSteam ('h_ps', p9, s9_id);
Wtur2_id = m_w*(h8-h9_id);

n_tur = 0.45; %***** variar con T********
Wtur2_real = Wtur2_id*n_tur; 
Wtur2_real_hp = Wtur2_real*1000*w_hp;

h9 = h8-(Wtur2_real/m_w);
s9 = XSteam ('s_ph', p9, h9);
T9 = XSteam ('T_ps', p9, s9);
x9 = XSteam ('x_ph', p4, h9);

incr_p2 = p8-p9;



%% RENDIMIENTO 

n_ciclo_alta = Wtur1_real/(Qint1*10^-3+Wb1_real);
n_ciclo_baja = Wtur2_real/(Qint2*10^-3+Qint3*10^-3+Wb2_real);

n_ciclo_total = (Wtur1_real+Wtur2_real)/(Qint3*10^-3+Wb2_real+Qint1*10^-3+Wb1_real)

%% POTENCIA TOTAL

WT_hp = Wtur1_real_hp + Wtur2_real_hp - Wb1_real_hp - Wb2_real_hp

%% DIAGRAMA TS

%% Curva liquido saturado 

TL_sat = linspace ( 0, 500);
SL_sat = zeros (1,100);

for nn = 1:100
    
  SL_sat(1,nn) =   XSteam ('sL_T',TL_sat(nn));
    
end

plot (SL_sat,TL_sat,'k')
hold on


%% Curva vapor saturado

TV_sat = linspace ( 0, 500);
SV_sat = zeros (1,100);

for nn = 1:100
    
  SV_sat(1,nn) =   XSteam ('sV_T',TV_sat(nn));
    
end

plot (SV_sat,TV_sat,'k')
hold on

%% Ciclo de alta temp.

plot (s1, T1,'r*')
hold on
plot (s2, T2,'r*')
hold on
plot (s3, T3,'r*')
hold on
plot (s4, T4,'r*')
hold on
plot (s4b, T4b, 'r*')
hold on


s1_2 = linspace (s1,s2);
h1_2 = linspace (h1,h2);
Y1_2 = zeros (1,100);

for nn = 1:100
    
  Y1_2(1,nn) =   XSteam ('T_hs',h1_2(nn), s1_2(nn));
    
end

plot (s1_2,Y1_2, 'c') 
hold on


s2_3 = linspace (s2,s3);
p2_3 = linspace (p2,p3);
Y2_3 = zeros (1,100);

for nn = 1:100
    
  Y2_3(1,nn) =   XSteam ('T_ps',p2_3(nn), s2_3(nn));
    
end

plot (s2_3,Y2_3, 'c') 
hold on


s3_4 = linspace (s3,s4);
h3_4 = linspace (h3,h4);
Y3_4 = zeros (1,100);

for nn = 1:100
    
  Y3_4(1,nn) =   XSteam ('T_hs',h3_4(nn), s3_4(nn));
    
end

plot (s3_4,Y3_4, 'c') 
hold on

s4_4b = linspace (s4,s4b);
p4_4b = linspace (p4,p4b);
Y4_4b = zeros (1,100);

for nn = 1:100
    
  Y4_4b(1,nn) =   XSteam ('T_ps',p4_4b(nn), s4_4b(nn));
    
end

plot (s4_4b,Y4_4b, 'c') 
hold on

s4b_1 = linspace (s4b,s1);
p4b_1 = linspace (p4b,p1);
Y4b_1 = zeros (1,100);

for nn = 1:100
    
  Y4b_1(1,nn) =   XSteam ('T_ps',p4b_1(nn), s4b_1(nn));
    
end

plot (s4b_1,Y4b_1, 'c') 
hold on

%% ciclo baja temperatura

plot (s5, T5,'m*')
hold on
plot (s6, T6,'m*')
hold on
plot (s7, T7,'m*')
hold on
plot (s8, T8,'m*')
hold on
plot (s9, T9,'m*')
hold on


s5_6 = linspace (s5,s6);
h5_6 = linspace (h5,h6);
Y5_6 = zeros (1,100);

for nn = 1:100
    
  Y5_6(1,nn) =   XSteam ('T_hs',h5_6(nn), s5_6(nn));
    
end

plot (s5_6,Y5_6, 'b') 
hold on


s6_7 = linspace (s6,s7);
p6_7 = linspace (p6,p7);
Y6_7 = zeros (1,100);

for nn = 1:100
    
  Y6_7(1,nn) =   XSteam ('T_ps',p6_7(nn), s6_7(nn));
    
end

plot (s6_7,Y6_7, 'b') 
hold on



s7_8 = linspace (s7,s8);
p7_8 = linspace (p7,p8);
Y7_8 = zeros (1,100);

for nn = 1:100
    
  Y7_8(1,nn) =   XSteam ('T_ps',p7_8(nn), s7_8(nn));
    
end

plot (s7_8,Y7_8, 'b') 
hold on


s8_9 = linspace (s8,s9);
h8_9 = linspace (h8,h9);
Y8_9 = zeros (1,100);

for nn = 1:100
    
  Y8_9(1,nn) =   XSteam ('T_hs',h8_9(nn), s8_9(nn));
    
end

plot (s8_9,Y8_9, 'b') 
hold on


s9_5 = linspace (s9,s5);
p9_5 = linspace (p9,p5);
Y9_5 = zeros (1,100);

for nn = 1:100
    
  Y9_5(1,nn) =   XSteam ('T_ps',p9_5(nn), s9_5(nn));
    
end

plot (s9_5,Y9_5, 'b') 
hold on
