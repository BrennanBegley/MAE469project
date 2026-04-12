sunmucononical = 1;
SunTUtodays= 58.13; %Days
SunAU_TUtokm_s= 29.79; %km/s
AUtoKm= 1.495*10^8; %Km

J2000epoch=juliandate(datetime(2000,1,1,11,58,00));
dur.Format = 'd';
plotdate=juliandate(datetime(2026,12,25,23,12,00));
TOF=linspace(169.5419,180,1000)
Mercuryproblem4ar=zeros(3,length(TOF))
for i=1:length(TOF)
    [Mercuryproblem4ar(:,i),Mercuryproblem4av]=universalTOF(sunmucononical,TOF(i),mercuryj2000r,mercuryj2000v)
end
plot3(Mercuryproblem4ar(1,:),Mercuryproblem4ar(2,:),Mercuryproblem4ar(3,:))