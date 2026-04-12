%%% ploting function

<<<<<<< Updated upstream
function planet_plot(Mercury_r,Venus_r, Earth_r, Mars_r, Jupiter_r, Saturn_r, Uranus_r,  Neptune_r,  Pluto_r)% Spacecraft_x, Spacescraft_y, Spacecraft_z)
    clf
    
    
=======
function planet_plot(Mercury_r,Mercury_v, Venus_r, Earth_r, Mars_r, Jupiter_r, Saturn_r, Uranus_r,  Neptune_r,  Pluto_r, Planet_TOF_inital, Planet_TOF_final, Planet_TOF_res )% Spacecraft_x, Spacescraft_y, Spacecraft_z)
mu = 1;
Mercury = struct( ...
    'a', 0.387099, ...
    'e', 0.205631, ...
    'inc', deg2rad(7.00487), ...
    'OMEGA', deg2rad(48.33167), ...
    'omega', deg2rad(29.12478), ...
    'theta', deg2rad(174.7944), ...
    'mu', 22031.868551); % km^3/s^2
    clf
    [mercuryj2000r,mercuryj2000v]=posandvelvector(Mercury,mu);
    sunmucononical = 1;
SunTUtodays= 58.13; %Days
SunAU_TUtokm_s= 29.79; %km/s
AUtoKm= 1.495*10^8; %Km

J2000epoch=juliandate(datetime(2000,1,1,11,58,00));
dur.Format = 'd';
plotdate=juliandate(datetime(2026,12,25,23,12,00));
TOF=linspace(169.5419,180,1000)
Mercury_orbit_r=zeros(3,length(TOF));
Mercury_orbit_v=zeros(3,length(TOF));
for i=1:length(TOF)
    [Mercury_orbit_r(:,i),Mercury_orbit_v]=universalTOF(sunmucononical,TOF(i),Mercury_r,Mercury_v);
end

>>>>>>> Stashed changes
    
    plot3(0,0,0,'oy',DisplayName="Sun") %sun
    hold on
    %plot locations of all the planets
    plot3(Mercury_r(1),Mercury_r(2),Mercury_r(3),"o" ,DisplayName="Mercury")
    plot3(Venus_r(1),Venus_r(2),Venus_r(3),"o" ,DisplayName="Venus")
    plot3(Earth_r(1),Earth_r(2),Earth_r(3),"o" ,DisplayName="Earth")
    plot3(Mars_r(1),Mars_r(2),Mars_r(3),"o" ,DisplayName="Mars")
    plot3(Jupiter_r(1), Jupiter_r(2), Jupiter_r(3),"o" ,DisplayName="Jupiter")
    plot3(Saturn_r(1),Saturn_r(2),Saturn_r(3),"o" ,DisplayName="Saturn")
    plot3(Uranus_r(1),Uranus_r(2),Uranus_r(3),"o" ,DisplayName="Uranus")
    plot3(Neptune_r(1),Neptune_r(2),Neptune_r(3),"o" ,DisplayName="Neptune")
<<<<<<< Updated upstream
    plot3( Pluto_r(1), Pluto_r(2), Pluto_r(3),"o" ,DisplayName=" Pluto")
=======
    plot3(Pluto_r(1),Pluto_r(2),Pluto_r(3),"o" ,DisplayName="pluto")
    
    %plot orbits
    plot3(Mercury_orbit_r(1,:),Mercury_orbit_r(2,:),Mercury_orbit_r(3,:))
>>>>>>> Stashed changes
    hold off
    legend
