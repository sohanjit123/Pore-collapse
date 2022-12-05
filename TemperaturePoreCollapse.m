clc,clear all

%%%Dynamic pore collapse
d=1e-9; %nanosecond
alpha0=1.3;
a0=20e-6; %initial size

Y1=400e6; %Yield 
rho=8.93e3; % density
Tm=1356;

 B=5765.26;
 eta_m=2e-3;
 C=385;

T0=293;

 Pcrit=(2*Y1/3)*log(alpha0/(alpha0-1)); 
 Prate=0.1 * 1e8; 

 Pload=Prate/1000;
 P=0:Pload:1.5*Pcrit;
% 
%   yield(T0,Y1,Tm)
%   visc(T0,Tm,B,eta_m)

 alpha=zeros(length(P),1);
 time=zeros(length(P),1);
 a=zeros(length(P),1);
 b=zeros(length(P),1);
  T=zeros(length(P),1);


 alpha(1)=alpha0;
 alpha(2)=alpha(1);
 time(1)=0;
 T(1)=T0;
a(1)=a0;
a(2)=a(1);
b(1)=bee(a(1),a0,alpha(1));
r0=(a(1)+b(1))/2;

for k=2:length(P)


         %Temperature loop%
 
   time(k)=k/1000;
  adot=(a(k)-a(k-1))/d;
   r=arr(r0,a0,a(k));
   rdotbyr=a(k)*a(k)*(adot)/(r^3);

   const=(d/rho/C);

   first_term=-2*yield(T(k-1),Y1,Tm)*rdotbyr;
   second_term=12*visc(T(k-1),Tm,B,eta_m)*rdotbyr*rdotbyr;

   T(k)=T(k-1)+const*(first_term+second_term);
   
   

%    alpha(k)=1+((alpha0-1)*(a(k)^3)/(a0^3));
%size loop
     
     b(k)=bee(a(k),a0,alpha(k));
   term1=2*yield(T(k),Y1,Tm)*log(b(k)/a(k));
   term2=-6*visc(T(k),Tm,B,eta_m)*a(k)*a(k)*adot*((1/(a(k)^3))-(1/(b(k)^3)));
   term3=rho*(adot^2)/2*(1-((a(k)^4)/(b(k)^4)));

    term4=term1+term2+term3-P(k);
    term4=term4/rho/((a(k)/b(k))-1);

    term4=term4-2*adot*adot/a(k);
    
    a(k+1)=term4*d*d+2*a(k)-a(k-1);

    alpha(k+1)= 1+((alpha0-1)*(a(k+1)^3)/(a0^3));

    if alpha(k)<1
        k
        break;
    end


end

% plot(time,alpha(1:end-1),'LineWidth',3)
% xlabel('Time (nanoseconds)')
% ylabel('Distension \alpha')
% ax=gca;
% ax.FontSize = 16.0;
% ax.LineWidth = 1.0;


time=time(1:k);
alpha=alpha(1:k);
plot(time,alpha,'LineWidth',3)

title('Dynamic Pore Collapse models')

xlabel('Time (nanoseconds)')
ylabel('Distension \alpha')
ax=gca;
ax.FontSize = 16.0;
ax.LineWidth = 1.0;

title('Dynamic Pore Collapse models')



function r=arr(r0,a0,a)

r=(r0^3-a0^3+a^3)^(1/3);

end

function b = bee(a,a0,alpha0)
b = (a^3+(a0^3/(alpha0-1)))^(1/3);
end

function Y=yield(T,Y1,Tm)
Y=Y1*(1-T/Tm);
end

function eta=visc(T,Tm,B,eta_m)
eta= eta_m*exp(B*((1/T)-(1/Tm)));
end