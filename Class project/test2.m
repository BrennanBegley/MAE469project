sunmucononical = 1;
SunTUtodays= 58.13; %Days
SunAU_TUtokm_s= 29.79; %km/s
AUtoKm= 1.495*10^8; %Km

J2000epoch=juliandate(datetime(2000,1,1,11,58,00));
dur.Format = 'd';
plotdate=juliandate(datetime(2026,12,25,23,12,00));
TOF=linspace(169.5419,180,100)

[Mercuryproblem4ar,Mercuryproblem4av]=universalTOF(sunmucononical,TOF,mercuryj2000r,mercuryj2000v);
