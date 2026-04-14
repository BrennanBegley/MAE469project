%%
% note depart date needs to be passed in as a julian date 
function [earthdepartrijk,earthdepartvijk,marsdepartrijk,marsdepartvijk,marsarrivalrijk,marsarivalvijk,transferorbitearthtomarsstruc,earthtotransferorbit1deltaV,delmarsflyby,marsdeparthyperbolicexcessvelocity,marsarivalhyperbolicexcessvelocity,marsflybydepartscvelocityijk,transferscv1ijk,transferscv2ijk]=eathdepatruetomars(SUNAUtokm,Sunstruc,shortlong,earthj2000rijk,earthj2000vijk,marsj2000rijk,marsj2000vijk,departdate,transferorbit1TOF,earthstuc,marsstruc,earthcircularorbitalt,sunmucononical,AUTUtoKms,j2000date,hyperbloicradius )
    %On the departure date 
    departTOF=(departdate-j2000date)/58.13;
    [earthdepartrijk,earthdepartvijk]=universalTOF(sunmucononical,departTOF,earthj2000rijk,earthj2000vijk); %AU AU/TU
    [marsdepartrijk,marsdepartvijk]=universalTOF(sunmucononical,departTOF,marsj2000rijk,marsj2000vijk); %AU AU/TU

    %mars position on arival date determined by TOF
    [marsarrivalrijk,marsarivalvijk]=universalTOF(sunmucononical,transferorbit1TOF,marsdepartrijk,marsdepartvijk); %AU AU/TU

    %Guass and SC needed velocities at departure and arrival
    [transferscv1ijk,transferscv2ijk,transferorbitearthtomarsstruc,transferorbitscr1ijk,transferorbitscr2ijk]= gaussspeedsandtransferorbitorbitalelements(earthdepartrijk,marsarrivalrijk,transferorbit1TOF,sunmucononical,shortlong);

    %calculate the earth hyperbolic excess velocity
    earthtotransferorbit1hyperbolicexess=(transferscv1ijk-earthdepartvijk)*AUTUtoKms;

    % calculate the earth needed delta v for circular orbit to heliocentric
    [earthtotransferorbit1deltaV]= deltavforcirculartohyperbolic(earthcircularorbitalt,earthstuc,earthtotransferorbit1hyperbolicexess);

    % Calculate the hyperbolic turn angle for mars flyby
    [delmarsflyby,marsdeparthyperbolicexcessvelocity,marsarivalhyperbolicexcessvelocity,marsflybydepartscvelocityijk]=hyperbolicturnangle(transferscv2ijk,marsarivalvijk,AUTUtoKms,hyperbloicradius,marsstruc,Sunstruc,SUNAUtokm);

end
