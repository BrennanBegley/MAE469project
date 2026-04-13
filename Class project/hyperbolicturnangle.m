%%
%the hyperbloicradius is the alt above the planet surface
function [del,vdscinf,vascinf,vscdeparthelocentric]=hyperbolicturnangle(scarivalvijk,planetarivalvijk,AUTUtoKms,hyperbloicradius,planetstruc)
    
    % hyperbolic excess speed at planet 
    scarivalv=scarivalvijk*AUTUtoKms;         %km/s spacecraft velocity at arrival
    planetarivalv=planetarivalvijk*AUTUtoKms; %km/s planet velocity at arrival
    vascinf=scarivalv-planetarivalv; %km/s
    vascinfnorm=norm(vascinf);

    % print these values
    fprintf('scarivalv (km/s): [%g %g %g]\n', scarivalv);
    fprintf('planetarivalv (km/s): [%g %g %g]\n', planetarivalv);
    fprintf('vascinf (km/s): [%g %g %g]\n', vascinf);
    fprintf('vascinfnorm (km/s): %g\n', vascinfnorm);

    %hyperbolic turn angle
    disp(hyperbloicradius+planetstruc.r)
    del= 2*asin(1/(1+((hyperbloicradius+planetstruc.r)*vascinfnorm^2)/planetstruc.mu))

    %hyperbolic transformation matrix
    H=[cos((del)), -sin((del)), 0; sin((del)), cos((del)), 0; 0,0,1];
    
    %hyperbolic increase speed
    vdscinf=H*vascinf;
    
    % velocity of the spacecraft in heliocentric reference 
    vscdeparthelocentric=(vdscinf+planetarivalvijk)/AUTUtoKms;
end