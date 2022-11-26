%Simplified Spherical model%
clc,clear all

%%%Simplified Static pore collapse
d=1e-9; %microsecond
alpha0=1.1;
 a0=20e-6; %initial size


 Y=276e6; %Yield 
 rho=2.7e3; % density

  Pcrit=(2*Y/3)*log(alpha0/(alpha0-1)); 
  Prate=0.2 * 1e8; 
  P=0:Prate/100:3*Pcrit;
  alpha=zeros(length(P),1);

  for i=1:length(alpha)

    if P(i)<Pcrit

        alpha(i)=alpha0;
    else

        alpha(i)=1/(1-exp(-3*P(i)/2/Y));
    end
  end

plot(alpha,'LineWidth',3)
ax = gca; 
ax.FontSize = 16; 
xlabel('timesteps','FontSize',16,'FontWeight','bold')
ylabel('\alpha','FontSize',16,'FontWeight','bold')


%Simplified dynamic pore collapse
% 
% 
 tau=sqrt(rho*a0*a0/(3*Y*((alpha0-1)^(2/3)))); 
 alpha_d=zeros(length(P),1);
 count=0;
 for i=1:length(alpha_d)
 
         
 
     if P(i)<=Pcrit
         alpha_d(i)=alpha0;
        count=count+1;
    end
end

for k = count:(length(P)-1)
 NonD=d*d/tau/tau/Y;
pterm=P(k)-(2*Y/3)*log(alpha_d(k)/(alpha_d(k)-1));
f1=tara(alpha_d(k),1);
f4=tara(alpha_d(k),4);
dif=alpha_d(k)-alpha_d(k-1);
term1= (1/6)*dif*dif*f4;
term2=NonD*pterm;

term3=(1/f1)*(term1-term2);

alpha_d(k+1)=2*alpha_d(k)-alpha_d(k-1)+term3;

 end
% 
hold on
plot(alpha_d,'LineWidth',3)

legend('Static','dynamic')
title('Simplified Pore Collapse models')


function y = tara(a,n)

y= (a-1)^(-n/3)-a^(-n/3);

end


% 
% 
% 
