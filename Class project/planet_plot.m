function planet_plot(Mercury_r, Mercury_v, Venus_r, Venus_v, Earth_r, Earth_v, ...
    Mars_r, Mars_v, Jupiter_r, Jupiter_v, Saturn_r, Saturn_v, Uranus_r, Uranus_v, ...
    Neptune_r, Neptune_v, Pluto_r, Pluto_v, Planet_TOF_initial, Planet_TOF_final, Planet_TOF_res)
%PLANET_PLOT Plot the heliocentric orbits of all planets (plus Pluto).
%   The function propagates each planet over the requested time span using
%   universalTOF, then plots the current position and traced orbit in 3D.

    clf
    muSunCanonical = 1;
    TOF = linspace(Planet_TOF_initial, Planet_TOF_final, Planet_TOF_res);

    % Preallocate propagated position histories.
    Mercury_orbit_r = zeros(3, length(TOF)); Mercury_orbit_v = zeros(3, length(TOF));
    Venus_orbit_r   = zeros(3, length(TOF)); Venus_orbit_v   = zeros(3, length(TOF));
    Earth_orbit_r   = zeros(3, length(TOF)); Earth_orbit_v   = zeros(3, length(TOF));
    Mars_orbit_r    = zeros(3, length(TOF)); Mars_orbit_v    = zeros(3, length(TOF));
    Jupiter_orbit_r = zeros(3, length(TOF)); Jupiter_orbit_v = zeros(3, length(TOF));
    Saturn_orbit_r  = zeros(3, length(TOF)); Saturn_orbit_v  = zeros(3, length(TOF));
    Uranus_orbit_r  = zeros(3, length(TOF)); Uranus_orbit_v  = zeros(3, length(TOF));
    Neptune_orbit_r = zeros(3, length(TOF)); Neptune_orbit_v = zeros(3, length(TOF));
    Pluto_orbit_r   = zeros(3, length(TOF)); Pluto_orbit_v   = zeros(3, length(TOF));

    for k = 1:length(TOF)
        [Mercury_orbit_r(:,k), Mercury_orbit_v(:,k)] = universalTOF(muSunCanonical, TOF(k), Mercury_r, Mercury_v);
        [Venus_orbit_r(:,k),   Venus_orbit_v(:,k)]   = universalTOF(muSunCanonical, TOF(k), Venus_r,   Venus_v);
        [Earth_orbit_r(:,k),   Earth_orbit_v(:,k)]   = universalTOF(muSunCanonical, TOF(k), Earth_r,   Earth_v);
        [Mars_orbit_r(:,k),    Mars_orbit_v(:,k)]    = universalTOF(muSunCanonical, TOF(k), Mars_r,    Mars_v);
        [Jupiter_orbit_r(:,k), Jupiter_orbit_v(:,k)] = universalTOF(muSunCanonical, TOF(k), Jupiter_r, Jupiter_v);
        [Saturn_orbit_r(:,k),  Saturn_orbit_v(:,k)]  = universalTOF(muSunCanonical, TOF(k), Saturn_r,  Saturn_v);
        [Uranus_orbit_r(:,k),  Uranus_orbit_v(:,k)]  = universalTOF(muSunCanonical, TOF(k), Uranus_r,  Uranus_v);
        [Neptune_orbit_r(:,k), Neptune_orbit_v(:,k)] = universalTOF(muSunCanonical, TOF(k), Neptune_r, Neptune_v);
        [Pluto_orbit_r(:,k),   Pluto_orbit_v(:,k)]   = universalTOF(muSunCanonical, TOF(k), Pluto_r,   Pluto_v);
    end

    plot3(0, 0, 0, 'oy', DisplayName='Sun')
    hold on

    % Plot current planet locations.
    plot3(Mercury_r(1), Mercury_r(2), Mercury_r(3), 'o', DisplayName='Mercury')
    plot3(Venus_r(1),   Venus_r(2),   Venus_r(3),   'o', DisplayName='Venus')
    plot3(Earth_r(1),   Earth_r(2),   Earth_r(3),   'o', DisplayName='Earth')
    plot3(Mars_r(1),    Mars_r(2),    Mars_r(3),    'o', DisplayName='Mars')
    plot3(Jupiter_r(1), Jupiter_r(2), Jupiter_r(3), 'o', DisplayName='Jupiter')
    plot3(Saturn_r(1),  Saturn_r(2),  Saturn_r(3),  'o', DisplayName='Saturn')
    plot3(Uranus_r(1),  Uranus_r(2),  Uranus_r(3),  'o', DisplayName='Uranus')
    plot3(Neptune_r(1), Neptune_r(2), Neptune_r(3), 'o', DisplayName='Neptune')
    plot3(Pluto_r(1),   Pluto_r(2),   Pluto_r(3),   'o', DisplayName='Pluto')

    % Plot propagated orbital traces.
    plot3(Mercury_orbit_r(1,:), Mercury_orbit_r(2,:), Mercury_orbit_r(3,:), '-w', HandleVisibility='off')
    plot3(Venus_orbit_r(1,:),   Venus_orbit_r(2,:),   Venus_orbit_r(3,:),   '-w', HandleVisibility='off')
    plot3(Earth_orbit_r(1,:),   Earth_orbit_r(2,:),   Earth_orbit_r(3,:),   '-w', HandleVisibility='off')
    plot3(Mars_orbit_r(1,:),    Mars_orbit_r(2,:),    Mars_orbit_r(3,:),    '-w', HandleVisibility='off')
    plot3(Jupiter_orbit_r(1,:), Jupiter_orbit_r(2,:), Jupiter_orbit_r(3,:), '-w', HandleVisibility='off')
    plot3(Saturn_orbit_r(1,:),  Saturn_orbit_r(2,:),  Saturn_orbit_r(3,:),  '-w', HandleVisibility='off')
    plot3(Uranus_orbit_r(1,:),  Uranus_orbit_r(2,:),  Uranus_orbit_r(3,:),  '-w', HandleVisibility='off')
    plot3(Neptune_orbit_r(1,:), Neptune_orbit_r(2,:), Neptune_orbit_r(3,:), '-w', HandleVisibility='off')
    plot3(Pluto_orbit_r(1,:),   Pluto_orbit_r(2,:),   Pluto_orbit_r(3,:),   '-w', HandleVisibility='off')

    legend
    title('Solar System on December 25, 2026 at 11:12 pm UTC')
    xlabel('x (AU)'); ylabel('y (AU)'); zlabel('z (AU)')
    grid on
    axis equal
    view(3)
    hold off
end
