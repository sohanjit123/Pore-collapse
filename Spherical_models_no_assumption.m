clc,clear all

%%%Static pore collapse
d=1; %nanosecond
alpha0=1.5;
 a0=20e-6; %initial porosity

G=25e9; 
Y=276e6; %Yield 
rho=2.7e3; % density

  Pcrit=(2*Y/3)*log(alpha0/(alpha0-1)); 
  Prate=0.2*1e8; 



 alpha1=(2*G*alpha0+Y)/(2*G+Y);
 P1 = 2*Y/3/alpha1;

 alpha2=(2*G*alpha0)/(2*G+Y);
 P2=2/3*Y*log(alpha2/(alpha2-1));

  P=0:Prate/100:2*P2;
    alpha=zeros(length(P),1);

%  
for i=1:length(P)

    i

    if P(i)<=P1

       

%         eqn = 4*G*(alpha0-alp)/(3*alp*(alp-1))==P(i);

%         alph=solve(eqn,alp);
%         alpha(i)=double(max(alph));

     alpha(i)=fzero(@(alp) 4*G*(alpha0-alp)/(3*alp*(alp-1))-P(i),alpha0 );

      elseif P(i)>P1 && P(i)<=P2
% %            
%          term1=2*G*(alpha0-alp)/(Y*alp);
%          term2=2*G*(alpha0-alp)/(Y*alp-Y);
% 
%          eqn = 2/3*Y*(1-term1+log(term2))==P(i);
%            alph=vpasolve(eqn,alp);
%         alpha(i)=double(max(alph));

   alpha(i)=fzero(@(alpy) 2/3*Y*(1-(2*G*(alpha0-alpy)/(Y*alpy))+log(2*G*(alpha0-alpy)/(Y*alpy-Y)))-P(i),alpha1);

    else
        
        alpha(i)=1/(1-exp(-3*P(i)/2/Y));
    end

end

plot(alpha)
plot(alpha,'LineWidth',3)
ax = gca; 
ax.FontSize = 16; 
xlabel('timesteps','FontSize',16,'FontWeight','bold')
ylabel('\alpha','FontSize',16,'FontWeight','bold')