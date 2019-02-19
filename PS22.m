%% Problem set 2
close all;
clear all;
clc;

%Global variables
global w1 w2 w3 rx1 rx2 rx3 Kl1 Kl2 Kl3 tl1 tl2 tl3 Sxp Rlt Dx Dl mu

%Basic Parameters
Dt = 30; %Doubling time (min)
DW = 0.3; %Percentage of dry mass per cell
Gc = 200; %Copies per cell
mRNA_h = 300/60; %mRNA half-life (min)
prot_h = 70*60; %Protein half-life (min)
Cv = 9*10^(-17); %Cell volume
Cm = 2.8*10^(-13); %Cell mass
Cc = 5*10^7; %Cell concentration (cell/mL)
Kep = 49*60; %Elongation rate (nts/min)
Klp = 16.5*60; %Translation rate (aa/min)
Ki = 0.024*60; %Initiation rate (1/min)
Rxt = 1.083*10^(-6); %Total RNAP concentration (M)
Lx1 = 1200; %gene1 length (nts)
Lx2 = 2400; %gene2 length (nts)
Lx3 = 600; %gene3 length (nts)
Ll1 = 400; %Protein 1 length (AA)
Ll2 = 800; %Protein 2 length (AA)
Ll3 = 200; %Protein 3 length (AA)
w1 = 0.26;
w2 = 30;
w3 = 30;
Av = 6.023*10^23; %Avogadro number
Rib = 20100; %Number of ribosomes per cell

%Compound Parameters
Rxt = Rxt*Cv*DW/Cm; %RNAP concentration (mol/gDW)
Rlt = (Rib*Cc*1000/Av)*(Cv*DW/Cm); %Ribosomes concentration (mol/gDW)
Dx = log(2)/mRNA_h; %Degradation rate of mRNA
Dl = log(2)/prot_h; %Degradation rate of proteins
mu = log(2)/Dt; %Dilution factor

Sxp = (1.04*Ki*Cv*DW/Cm)*10^(-6); %Saturation constant (mol/gDW)
Ke1 = Kep/Lx1; %Elongation rate for mRNA1
Ke2 = Kep/Lx2; %Elongation rate for mRNA2
Ke3 = Kep/Lx3; %Elongation rate for mRNA3
tx1 = Ke1/Ki; %Tau for mRNA1
tx2 = Ke2/Ki; %Tau for mRNA2
tx3 = Ke3/Ki; %Tau for mRNA3
MW1 = Lx1*607.4+157.9; %Molecular Weight of gene 1
MW2 = Lx2*607.4+157.9; %Molecular Weight of gene 2
MW3 = Lx3*607.4+157.9; %Molecular Weight of gene 3
Gp1 = (Cc*1000*Gc*50*10^(-9)/MW1)*(Cv*DW/Cm); %Gp for gene 1
Gp2 = (Cc*1000*Gc*50*10^(-9)/MW2)*(Cv*DW/Cm); %Gp for gene 2
Gp3 = (Cc*1000*Gc*50*10^(-9)/MW3)*(Cv*DW/Cm); %Gp for gene 3

rx1 = Ke1*Rxt*(Gp1/(Sxp*tx1+Gp1*tx1+Gp1)); %Transcription rate for gene 1
rx2 = Ke2*Rxt*(Gp2/(Sxp*tx2+Gp2*tx2+Gp2)); %Transcription rate for gene 2
rx3 = Ke3*Rxt*(Gp3/(Sxp*tx3+Gp3*tx3+Gp3)); %Transcription rate for gene 3

Kl1 = Klp/Ll1; %Translation rate protein 1
Kl2 = Klp/Ll2; %Translation rate protein 2
Kl3 = Klp/Ll3; %Translation rate protein 3
tl1 = Kl1/Ki; %Tau for protein 1
tl2 = Kl2/Ki; %Tau for protein 2
tl3 = Kl3/Ki; %Tau for protein 3

%Initial conditions
t_i = 0;
t_f = 460; %Final time (min)
step = 0.01;
t_span = t_i:step:t_f; %Time vector (min)
[m,n] = size(t_span); %Size of time
I = zeros(n+1,1);
I(1:60,1) = 10*10^(-3); %Inducer initial concentration (mol/gDW)
x0 = [0;0;0;0;0;0]; %Initial conditions for x vector

%% Case 1 - Normal circuit

[t,X] = ode15s(@(t,x) sys(t,x,I),t_span,x0);
X = X.*(10^(9));

figure(1)
q = plot(t_span,X(:,1),t_span,X(:,2),t_span,X(:,3));
q(1).LineWidth = 1.2;
q(1).Color = 'black';
q(1).LineStyle = '-';
q(2).LineWidth = 1.2;
q(2).Color = [0.4 0.4 0.5];
q(2).LineStyle = '--';
q(3).LineWidth = 1.2;
q(3).Color = [0.6 0.6 0.6];
q(3).LineStyle = '-.';
xlabel('Time [min]','fontweight','bold')
ylabel('Protein concentration [umol/gDW]','fontweight','bold')
legend('p1','p2','p3')

%% Case 2 - Broken circuit

[t,Y] = ode15s(@(t,x) sys2(t,x,I),t_span,x0);
Y = Y.*(10^(9));

figure(2)
q = plot(t_span,Y(:,1),t_span,Y(:,2),t_span,Y(:,3));
q(1).LineWidth = 1.2;
q(1).Color = 'black';
q(1).LineStyle = '-';
q(2).LineWidth = 1.2;
q(2).Color = [0.4 0.4 0.5];
q(2).LineStyle = '--';
q(3).LineWidth = 1.2;
q(3).Color = [0.6 0.6 0.6];
q(3).LineStyle = '-.';
xlabel('Time [min]','fontweight','bold')
ylabel('Protein concentration [umol/gDW]','fontweight','bold')
legend('p1','p2','p3')


%% Case 3 - Approximate solution

step = 0.01; %(min)
t_span = t_i:step:t_f; %Time vector (min)
[m,l] = size(t_span); %Extract the size of time
x = zeros(6,l);
r = zeros(6,1);
A = [-Dx 0 0 0 0 0;...
    0 -(Dx+mu) 0 0 0 0;...
    0 0 -(Dx+mu) 0 0 0;...
    0 0 0 -Dl 0 0;...
    0 0 0 0 -(Dl+mu) 0;...
    0 0 0 0 0 -(Dl+mu)];
S = eye(6,6);

A1 = exp(A)*step;
S1 = inv(A)*(A1-eye(6,6))*S;

n = 15;
K = 0.3*10^(-9);

for i=1:l-1
    if i<6000
        In = 10*10^(-3);
    else
        In = 0;
    end
    fi = (In^n)/(In^n+K^n); %Inducer function
    fp1 = (x(1,i)^n)/(x(1,i)^n+K^n); %P1 function
    fp2 = (x(2,i)^n)/(x(2,i)^n+K^n); %P2 function
    fp3 = (x(3,i)^n)/(x(3,i)^n+K^n); %P3 function

    ui = (w1+w2*fi)/(1+w1+w2*fi); %Control inducer function
    up13 = (w1+w2*fp1+w3*fp3)/(1+w1+w2*fp1+w3*fp3); %Control p2 function
    up12 = (w1+w2*fp1+w3*fp2)/(1+w1+w2*fp1+w3*fp2); %Control p3 function
%     up12 = (w1+w2*fp1)/(1+w1+w2*fp1); %Control p3 function broken
%     (uncomment to obtain broken graph)
    
    r(1,1) = rx1*ui;
    r(2,1) = rx2*up13;
    r(3,1) = rx3*up12;
    r(4,1) = Kl1*Rlt*(x(1,i)/(Sxp*tl1+x(1,i)*tl1+x(1,i))); 
    r(5,1) = Kl2*Rlt*(x(2,i)/(Sxp*tl2+x(2,i)*tl2+x(2,i))); 
    r(6,1) = Kl3*Rlt*(x(3,i)/(Sxp*tl3+x(3,i)*tl3+x(3,i))); 


    x(:,i+1) = A1*x(:,i)+S1*r(:,1);
end

x = x.*10^9;

figure(3)
q = plot(t_span,x(1,:),t_span,x(2,:),t_span, x(3,:));
q(1).LineWidth = 1.2;
q(1).Color = 'black';
q(1).LineStyle = '-';
q(2).LineWidth = 1.2;
q(2).Color = [0.4 0.4 0.5];
q(2).LineStyle = '--';
q(3).LineWidth = 1.2;
q(3).Color = [0.6 0.6 0.6];
q(3).LineStyle = '-.';
xlabel('Time [min]','fontweight','bold')
ylabel('Protein concentration [umol/gDW]','fontweight','bold')
legend('p1','p2','p3')
