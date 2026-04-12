
% expects angles in radians 
% R and V from Orbital elements 
function [rijk,vijk] = posandvelvector(planet,mu)
    a       = planet.a;
    e       = planet.e;
    inc     = planet.inc;
    OMEGA   = planet.OMEGA;
    omega   = planet.omega;
    trueanom   = planet.theta;
    %PQW plane 
    p=a*(1-e^2); 
    r=p/(1+e*cos(trueanom));
    rpqw= [r*cos(trueanom);r*sin(trueanom);0];
    vpqw= sqrt(mu/p)*[-sin(trueanom); e+cos(trueanom);0];
    %PQW to IJK matrix
    R = [cos(OMEGA)*cos(omega) - sin(OMEGA)*sin(omega)*cos(inc),  -cos(OMEGA)*sin(omega) - sin(OMEGA)*cos(omega)*cos(inc),  sin(OMEGA)*sin(inc);
         sin(OMEGA)*cos(omega) + cos(OMEGA)*sin(omega)*cos(inc),  -sin(OMEGA)*sin(omega) + cos(OMEGA)*cos(omega)*cos(inc), -cos(OMEGA)*sin(inc);
         sin(omega)*sin(inc),                                      cos(omega)*sin(inc),                                     cos(inc)           ];
    %IJK plane
    rijk=R*rpqw;
    vijk=R*vpqw;
end