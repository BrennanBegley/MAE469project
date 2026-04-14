function planet_plot(Mercury_r, Mercury_v, Venus_r, Venus_v, Earth_r, Earth_v, ...
    Mars_r, Mars_v, Jupiter_r, Jupiter_v, Saturn_r, Saturn_v, Uranus_r, Uranus_v, ...
    Neptune_r, Neptune_v, Pluto_r, Pluto_v, Planet_TOF_initial, Planet_TOF_final, Planet_TOF_res)
%PLANET_PLOT Plot all planets in dark mode with exactly one full orbit each.

    clf
    muSunCanonical = 1;

    % =========================
    % Force dark mode
    % =========================
    set(gcf, 'Color', [0.05 0.05 0.05])

    ax = gca;
    ax.Color = [0.02 0.02 0.02];
    ax.XColor = [1 1 1];
    ax.YColor = [1 1 1];
    ax.ZColor = [1 1 1];
    ax.GridColor = [0.6 0.6 0.6];
    ax.MinorGridColor = [0.4 0.4 0.4];

    hold on

    % =========================
    % Sun
    % =========================
    plot3(0, 0, 0, 'yo', ...
        'MarkerFaceColor', 'y', ...
        'MarkerSize', 10, ...
        'DisplayName', 'Sun')

    % =========================
    % Plot planets
    % =========================
    plotSinglePlanetOneOrbit(Mercury_r, Mercury_v, Planet_TOF_res, muSunCanonical, 'Mercury');
    plotSinglePlanetOneOrbit(Venus_r,   Venus_v,   Planet_TOF_res, muSunCanonical, 'Venus');
    plotSinglePlanetOneOrbit(Earth_r,   Earth_v,   Planet_TOF_res, muSunCanonical, 'Earth');
    plotSinglePlanetOneOrbit(Mars_r,    Mars_v,    Planet_TOF_res, muSunCanonical, 'Mars');
    plotSinglePlanetOneOrbit(Jupiter_r, Jupiter_v, Planet_TOF_res, muSunCanonical, 'Jupiter');
    plotSinglePlanetOneOrbit(Saturn_r,  Saturn_v,  Planet_TOF_res, muSunCanonical, 'Saturn');
    plotSinglePlanetOneOrbit(Uranus_r,  Uranus_v,  Planet_TOF_res, muSunCanonical, 'Uranus');
    plotSinglePlanetOneOrbit(Neptune_r, Neptune_v, Planet_TOF_res, muSunCanonical, 'Neptune');
    plotSinglePlanetOneOrbit(Pluto_r,   Pluto_v,   Planet_TOF_res, muSunCanonical, 'Pluto');

    % =========================
    % Perfect equal axes
    % =========================
    axis equal
    daspect([1 1 1])
    pbaspect([1 1 1])


    grid on
    view(3)

    xlabel('x (AU)', 'Color', 'w')
    ylabel('y (AU)', 'Color', 'w')
    zlabel('z (AU)', 'Color', 'w')
    title('Solar System on December 25, 2026 at 11:12 pm UTC', ...
        'Color', 'w')

    legend('TextColor', 'w', ...
           'Color', [0.1 0.1 0.1], ...
           'Location', 'bestoutside')

    hold off
end
