%%
%the hyperbloicradius is the alt above the planet surface
function [del,vdscinfIJK,vascinfijk,vscdeparthelocentric]=hyperbolicturnangle(scarivalvijk,planetarivalvijk,AUTUtoKms,hyperbloicradius,planetstruc,Sunstruc,SUNAUtokm)
    
    %calculate the body SOI
    PlanetSOI=norm(planetstruc.a)*(planetstruc.mass/Sunstruc.mass)^(2/5);

    %calculate the spacecraft normalized radial vector in the planet frame 
    SCmarsSOIframearivalrijk= -PlanetSOI*scarivalvijk/norm(scarivalvijk); %AU

    %calculate the SC hyperbolic excess velocity
    scarivalv=norm(scarivalvijk);         %AU/TU spacecraft velocity at arrival
    planetarivalv=norm(planetarivalvijk); %AU/TUplanet velocity at arrival
    vascinfijk=scarivalvijk-planetarivalvijk; %AU/TU
    
    %calculate the hyperbolic orbit OE 
    [a,e,enorm,i,RAAN,argumentperi,trueanom]=orbitalelementscalc(SCmarsSOIframearivalrijk*SUNAUtokm,vascinfijk*AUTUtoKms,planetstruc.mu);
    
    % Calculate the perafocal frame 
    % perafocal frame and spacecraft PQW velocity
    p_sc=a*(1-enorm^2); %AU
    fprintf('p_sc (AU): %g\n', p_sc/SUNAUtokm);
    
    r_scpqw=p_sc/(1+enorm*cos(trueanom)); %AU
    Rsc_pqwplanetframe=[r_scpqw*cos(trueanom); r_scpqw*sin(trueanom); 0];
    P_unitvector= [-sin(trueanom); 0; 0];
    Q_unitvector= [0; enorm+cos(trueanom); 0];
    V_scPQWplanetframe=sqrt(planetstruc.mu/(p_sc))*(P_unitvector+Q_unitvector); %km/s

    % Print header comment and the computed PQW vectors/values
    fprintf('%% Perifocal frame and spacecraft PQW velocity results\n');
    fprintf('Rsc_pqwplanetframe (km): [%g, %g, %g]\n', Rsc_pqwplanetframe);
    fprintf('V_scPQWplanetframe (km/s): [%g, %g, %g]\n', V_scPQWplanetframe);
    fprintf('p_sc (AU): %g, r_scpqw (AU): %g, eccentricity e: %g, true anomaly (rad): %g\n', p_sc, r_scpqw, enorm, trueanom);

    %calculate the hyperbolic turn 
    rp_sc_hyperbolic= planetstruc.r+hyperbloicradius; %km
    del= 2*asin(1/(1+((rp_sc_hyperbolic)*norm(vascinfijk*AUTUtoKms)^2)/planetstruc.mu)); %rads
    
    %hyperbolic transformation matrix
    H=[cos((del)), -sin((del)), 0; sin((del)), cos((del)), 0; 0,0,1];
    
    %hyperbolic speed PQW
    vdscinfpqw=H*V_scPQWplanetframe %km/s

    %PQW to IJK (angles are in radians)
    R1 = [cos(RAAN)*cos(argumentperi)-sin(RAAN)*sin(argumentperi)*cos(i), -cos(RAAN)*sin(argumentperi)-sin(RAAN)*cos(argumentperi)*cos(i),  sin(RAAN)*sin(i);
          sin(RAAN)*cos(argumentperi)+cos(RAAN)*sin(argumentperi)*cos(i), -sin(RAAN)*sin(argumentperi)+cos(RAAN)*cos(argumentperi)*cos(i), -cos(RAAN)*sin(i);
          sin(argumentperi)*sin(i),                                         cos(argumentperi)*sin(i),                                        cos(i)];
    
    %PQW to IJK 
    vdscinfIJK=R1*vdscinfpqw %km/s
    vscdeparthelocentric=(vdscinfIJK+planetarivalvijk*AUTUtoKms)/AUTUtoKms %AU/TU
end