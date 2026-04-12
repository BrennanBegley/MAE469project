%%% ploting function

function planet_plot(Mercury_r,Venus_r, Earth_r, Mars_r, Jupiter_r, Saturn_r, Uranus_r,  Neptune_r,  Pluto_r)% Spacecraft_x, Spacescraft_y, Spacecraft_z)
    clf
    
    
    
    plot3(0,0,0,'oy',DisplayName="Sun") %sun
    hold on

    plot3(Mercury_r(1),Mercury_r(2),Mercury_r(3),"o" ,DisplayName="Mercury")
    plot3(Venus_r(1),Venus_r(2),Venus_r(3),"o" ,DisplayName="Venus")
    plot3(Earth_r(1),Earth_r(2),Earth_r(3),"o" ,DisplayName="Earth")
    plot3(Mars_r(1),Mars_r(2),Mars_r(3),"o" ,DisplayName="Mars")
    plot3(Jupiter_r(1), Jupiter_r(2), Jupiter_r(3),"o" ,DisplayName="Jupiter")
    plot3(Saturn_r(1),Saturn_r(2),Saturn_r(3),"o" ,DisplayName="Saturn")
    plot3(Uranus_r(1),Uranus_r(2),Uranus_r(3),"o" ,DisplayName="Uranus")
    plot3(Neptune_r(1),Neptune_r(2),Neptune_r(3),"o" ,DisplayName="Neptune")
    plot3( Pluto_r(1), Pluto_r(2), Pluto_r(3),"o" ,DisplayName=" Pluto")
    hold off
    legend
