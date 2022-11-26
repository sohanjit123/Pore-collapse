clc,clear all

%%%Dynamic pore collapse
d=1e-9; %nanosecond
alpha0=1.1;
a0=20e-6; %initial size


Y=276e6; %Yield 
rho=2.7e3; % density

Pcrit=(2*Y/3)*log(alpha0/(alpha0-1)); 
Prate=0.2 * 1e8; 
P=0:Prate/1000:1.5*Pcrit;

 alpha_d=zeros(length(P),1);
G=25e9; 


alpha_d(1)=alpha0;
 alpha_d(2)=alpha_d(1);



 for k=2:length(P)


a=sqrt(a0*a0*(alpha_d(k)-1)/(alpha0-1));
b=sqrt(a0*a0*alpha_d(k)/(alpha0-1));
B=a0*a0*(alpha0-alpha_d(k))/(alpha0-1);

if (4*G*B/a/a)<=Y

    %Pure elastic%

rhs=P(k)-2*G*B*((1/a/a)-(1/b/b));
rhs=rhs/rho;
bdott=  -1*((alpha_d(k)-alpha_d(k-1))/d)*a0*a0/(alpha0-1);
deductterm=(bdott^2)/4*((1/b/b)-(1/a/a));
rhs=rhs-deductterm;
rhs=rhs*2/log(b/a);
alpha_d(k+1)=(rhs*d*d*(alpha0-1)/(-1*a0*a0))+2*alpha_d(k)-alpha_d(k-1);

elseif (4*G*B/a/a)>Y && (4*G*B/b/b)<Y

    %Elastic-plastic
c=sqrt(4*G*B/Y);

extraterm=Y*log(c/a);
rhs=P(k)-2*G*B*((1/c/c)-(1/b/b))-extraterm;
rhs=rhs/rho;
bdott=  -1*((alpha_d(k)-alpha_d(k-1))/d)*a0*a0/(alpha0-1);
deductterm=(bdott^2)/4*((1/b/b)-(1/a/a));
rhs=rhs-deductterm;
rhs=rhs*2/log(b/a);
alpha_d(k+1)=(rhs*d*d*(alpha0-1)/(-1*a0*a0))+2*alpha_d(k)-alpha_d(k-1);


else
    extraterm=Y*log(b/a);
rhs=P(k)-extraterm;
rhs=rhs/rho;
bdott=  -1*((alpha_d(k)-alpha_d(k-1))/d)*a0*a0/(alpha0-1);
deductterm=(bdott^2)/4*((1/b/b)-(1/a/a));
rhs=rhs-deductterm;
rhs=rhs*2/log(b/a);
alpha_d(k+1)=(rhs*d*d*(alpha0-1)/(-1*a0*a0))+2*alpha_d(k)-alpha_d(k-1);

    
    %fully plastic
end


 end



 hold on
plot(alpha_d,'LineWidth',3)


title('Dynamic Cylindrical Pore Collapse models')

