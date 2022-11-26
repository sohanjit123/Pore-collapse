%Simplified Spherical model%
clc,clear all

%%%Dynamic pore collapse
d=1e-9; %nanosecond
alpha0=1.5;
 a0=20e-6; %initial size


 Y=276e6; %Yield 
 rho=2.7e3; % density

  Pcrit=(2*Y/3)*log(alpha0/(alpha0-1)); 
  Prate=0.2 * 1e8; 
  P=0:Prate/1000:1.5*Pcrit;


    G=25e9; 

 alpha1=(2*G*alpha0+Y)/(2*G+Y);
 alpha2=(2*G*alpha0)/(2*G+Y);

 tau=sqrt(rho*a0*a0/(3*Y*((alpha0-1)^(2/3)))); 
 alpha_d=zeros(length(P),1);
 count=0;
    
 alpha_d(1)=alpha0;
 alpha_d(2)=alpha_d(1);


   NonD=d*d/tau/tau/Y;

 for k=2:length(alpha_d)
 
         
 
     if alpha_d(k)>alpha1


         sterm=P(k)-4*G*(alpha0-alpha_d(k))/(3*alpha_d(k)*(alpha_d(k)-1));
            f1=tara(alpha_d(k),1);
            f4=tara(alpha_d(k),4);
            dif=alpha_d(k)-alpha_d(k-1);
            term1= (1/6)*dif*dif*f4;
            term4=NonD*sterm;

            term3=(1/f1)*(term1-term4);

            alpha_d(k+1)=2*alpha_d(k)-alpha_d(k-1)+term3;
        
     elseif alpha_d(k)<=alpha1 && alpha_d(k)>alpha2
            
        sterm=P(k)-2/3*Y*(1-(2*G*(alpha0-alpha_d(k))/(Y*alpha_d(k)))+log(2*G*(alpha0-alpha_d(k))/(Y*alpha_d(k)-Y)));
          
          f1=tara(alpha_d(k),1);
            f4=tara(alpha_d(k),4);
            dif=alpha_d(k)-alpha_d(k-1);
            term1= (1/6)*dif*dif*f4;
            term4=NonD*sterm;

            term3=(1/f1)*(term1-term4);

            alpha_d(k+1)=2*alpha_d(k)-alpha_d(k-1)+term3;

     else

         
            
            sterm=P(k)-(2*Y/3)*log(alpha_d(k)/(alpha_d(k)-1));
            f1=tara(alpha_d(k),1);
            f4=tara(alpha_d(k),4);
             dif=alpha_d(k)-alpha_d(k-1);
            term1= (1/6)*dif*dif*f4;
            term2=NonD*sterm;

            term3=(1/f1)*(term1-term2);

            alpha_d(k+1)=2*alpha_d(k)-alpha_d(k-1)+term3;
       
    end
end

plot(alpha_d,'LineWidth',3)


title('Dynamic Pore Collapse models')


function y = tara(a,n)

y= (a-1)^(-n/3)-a^(-n/3);

end


