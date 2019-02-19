function f = sys(t,x,I)
global w1 w2 w3 rx1 rx2 rx3 Kl1 Kl2 Kl3 tl1 tl2 tl3 Sxp Rlt Dx Dl mu

%x = [m1;m2;m3;p1;p2;p3];
n = 25;
K = 0.3*10^(-9);

a = round(t);
In = I(a+1,1);
m1 = x(1);
m2 = x(2);
m3 = x(3);
p1 = x(4);
p2 = x(5);
p3 = x(6);

fi = (In^n)/(In^n+K^n); %Inducer function
fp1 = (m1^n)/(m1^n+K^n); %P1 function
fp2 = (m2^n)/(m2^n+K^n); %P2 function
fp3 = (m3^n)/(m3^n+K^n); %P3 function


ui = (w1+w2*fi)/(1+w1+w2*fi); %Control inducer function
up13 = (w1+w2*fp1+w3*fp3)/(1+w1+w2*fp1+w3*fp3); %Control p2 function
up12 = (w1+w2*fp1+w3*fp2)/(1+w1+w2*fp1+w3*fp2); %Control p3 function


rl1 = Kl1*Rlt*(m1/(Sxp*tl1+m1*tl1+m1))*10^3; %Translation rate for gene 1
rl2 = Kl2*Rlt*(m2/(Sxp*tl2+m2*tl2+m2))*10^3; %Translation rate for gene 2
rl3 = Kl3*Rlt*(m3/(Sxp*tl3+m3*tl3+m3))*10^3; %Translation rate for gene 3 


f = zeros(6,1);
f(1,1) = rx1*ui-m1*Dx;
f(2,1) = rx2*up13-m2*(Dx+mu);
f(3,1) = rx3*up12-m3*(Dx+mu);
f(4,1) = rl1-p1*Dl;
f(5,1) = rl2-p2*(Dl+mu);
f(6,1) = rl3-p3*(Dl+mu);

end