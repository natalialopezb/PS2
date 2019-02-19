%% Problem set 1.c
close all;
clear all;
clc;

%% Parameters
%Weight function
w1 = 0.26;
w2 = 300;
K = 0.3/1000; %In M
n = 1.5;

I_i = 0.0001/1000; %In M
I_f = 10/1000; %In M
I = I_i:0.001/1000:I_f;

f = (I.^(n))./(K^n+I.^n); %Inducer function

%Rates and initial values
%Basic parameters
L = 3075; %Gene length (nt)
Gc = 2500; %Copies per cell
ki = 0.024; %Initiation rate (s-1)
ka = 3.3*10^(-5); %Abortion rate (s-1)
Rx = 1.083*10^(-6); %Total RNAP per cell
deltaC = 5*10^7; %Cell concentration (cells/mL)
kes = 49; %Elongation rate (nts/s)

Par = [L Gc ki ka Rx deltaC kes]; %Vector with basic parameters
name = ["Gene  length", "Copies per cell", "Initiation rate", ...
    "Abortion rate", "Total RNAP", "Cell concentration","Elongation rate"];

%% Let's loop!

[a,stop] = size(Par);

for i=1:stop
    %Before changing things
    temp_Par = Par(1,i);
    ke = Par(1,7)/Par(1,1); %Specific elongation rate (s-1)
    Sx = 1.04*Par(1,3)*10^(-6); %Saturation constant (M)
    tx = (Par(1,4)*ke)/Par(1,3); %Dimentionless tau
    MWg = (Par(1,1)*607.4)+157.9; %Molecular Weight of plasmid (Da)
    Gp = (Par(1,6)*Par(1,2)*50*10^(-9))/MWg; %Gp in M. 50ng of plasmid are
                                      %suggested per transformation
    temp_rx = (ke*Par(1,5)*Gp)/(Sx*tx+Gp*tx+Gp);
    mp = temp_rx*((w1+w2.*f)./(1+w1+w2.*f)); %mRNA in M
    mp = mp.*10^6; %uM
    I = I.*10^3; %mM
    
    Par(1,i) = 10*Par(1,i); %Increase parameter 1 order of magnitude
    ke = 49/Par(1,1); %Specific elongation rate (s-1)
    Sx = 1.04*Par(1,3)*10^(-6); %Saturation constant (M)
    tx = (Par(1,4)*ke)/Par(1,3); %Dimentionless tau
    MWg = (Par(1,1)*607.4)+157.9; %Molecular Weight of plasmid (Da)
    Gp = (Par(1,6)*1000*Par(1,2)*50*10^(-9))/MWg; %Gp in M. 50ng of plasmid are
                                      %suggested per transformation
    in_rx = (ke*Par(1,5)*Gp)/(Sx*tx+Gp*tx+Gp);
    in_mp = in_rx*((w1+w2.*f)./(1+w1+w2.*f)); %mRNA in M
    in_mp = in_mp.*10^6; %uM
    
    Par(1,i) = Par(1,i)/100000; %Decrease parameter 1 order of magnitude
    ke = 49/Par(1,1); %Specific elongation rate (s-1)
    Sx = 1.04*Par(1,3)*10^(-6); %Saturation constant (M)
    tx = (Par(1,4)*ke)/Par(1,3); %Dimentionless tau
    MWg = (Par(1,1)*607.4)+157.9; %Molecular Weight of plasmid (Da)
    Gp = (Par(1,6)*1000*Par(1,2)*50*10^(-9))/MWg; %Gp in M. 50ng of plasmid are
                                      %suggested per transformation
    de_rx = (ke*Par(1,5)*Gp)/(Sx*tx+Gp*tx+Gp);
    de_mp = de_rx*((w1+w2.*f)./(1+w1+w2.*f)); %mRNA in M
    de_mp = de_mp.*10^6; %uM
    
    figure(i)
    p = semilogx(I,mp);
    hold on;
    q = semilogx(I,in_mp);
    r = semilogx(I,de_mp);
    p.Color = 'k';
    p.LineWidth = 0.85;
    q.Color = 'b';
    q.LineWidth = 0.85;
    r.Color = 'r';
    r.LineWidth = 0.85;
    xlabel('Inducer concentration [mM]','fontweight','bold');
    ylabel('mRNA concentration [uM]','fontweight','bold');
    conc1 = strcat('Increasing'," ",name(1,i));
    conc2 = strcat('Decreasing'," ",name(1,i));
    legend('Baseline',conc1,conc2);
    hold off;
    
    Par(1,i) = temp_Par;
end

%It looks like the [] part of r_x,p approaches 1 since Gp is much greater
%than Sx or tx, which means that the only parameters that affect this
%system are the elongation rate and RNAP concentration...