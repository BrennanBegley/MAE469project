%%
function [deltav]= deltavforcirculartohyperbolic(circularorbitalt,planetstruc,hyperbolicexcessvelocity)
    
    hyperbolicexcessvelocity=norm(hyperbolicexcessvelocity);
    
    %hyperbolic excess specific energy 
    hyperbolicorbitenegry=hyperbolicexcessvelocity^2/2; %km^2/s^2

    %hyperbolic excess velocity at circular PE
    vofhyperbolicatcircularpe= sqrt(2*(hyperbolicorbitenegry+planetstruc.mu/(planetstruc.r+circularorbitalt)));

    %circular orbit velocity 
    sccircorbitv=sqrt(planetstruc.mu/(planetstruc.r+circularorbitalt));

    %earth delta V
    deltav=vofhyperbolicatcircularpe-sccircorbitv;
end

