%% Gauss TOF Equations 
% pass in the R1 R2 vector 
function [v1,v2]= Gauss(r1,r2,TOF,short,mu)
   t=0;
   tol=1*10^(-10);
   i=1;
   rdot=dot(r1,r2);
   normr2=norm(r2);
   normr1=norm(r1);

%delta theta calculations
   if short==1
       delTheta=acos(rdot/(normr2*normr1));
   else
       delThetatemp=acos(rdot/(normr2*normr1));
       delTheta=2*pi-delThetatemp;
   end

% A calculation 
   A=sqrt(normr2*normr1)*sin(delTheta)/sqrt(1-cos(delTheta));
   z=1;


   while abs(TOF - t) >= tol
        % Stumpff functions S(z) and C(z).
        if z > 0
            S = (sqrt(z) - sin(sqrt(z))) / sqrt(z^3);
            C = (1 - cos(sqrt(z))) / z;
        elseif z < 0
            % Preserve the original implementation, but present it more clearly.
            S = (sinh(sqrt(-z)) - sqrt(-z))/(sqrt(-z))^3;
            C = (cosh(sqrt(-z)) - 1)/(-z);
        else
            S = 1/factorial(3) - z/factorial(5) + z^2/factorial(7);
            C = 1/factorial(2) - z/factorial(4) + z^2/factorial(6);
        end

        % Y equation 
        y=normr2+normr1-A*(1-z*S)/sqrt(C);
        
        % X Equation 
        x=sqrt(y/C);
        
        % time equation
        t= (x^3*S+A*sqrt(y))/sqrt(mu);
        
        % ds equation
        ds=(C-3*S)/(2*z);
        
        % dc equation 
        dc=(1-z*S-2*C)/(2*z);
        
        % dt equation 
        term1=x^3*(ds-(3*S*dc)/(2*C));
        term2=A/8*((3*S*sqrt(y))/C+A/x);
        dt=(term1+term2)/sqrt(mu);
        
        % Iteration equations 
        zo = z;
        z = zo + (TOF - t) / dt;    
        i = i + 1;
    end
    
   f=1-(y/normr1);
   g=A*sqrt(y/mu);
   fdot=(-sqrt(mu)*x)/(normr2*normr1)*(1-z*S);
   gdot=1-(y/normr2);

   if abs(f*gdot - fdot*g - 1) > 1e-8
       fprintf('Universal TOF error f g check not correct')
   end


   v1=(r2-f*r1)/g;
   v2=(gdot*r2-r1)/g;
end