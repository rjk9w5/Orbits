% Access Data sheet
function constants()
AU = 1.4960e8;

% sun constants
sun.axial_rotational_period_days = 27;
sun.radius_km = 696000;
sun.mu = 1.327e11;

save -mat-binary sun.mat sun

% moon constants
moon.axial_rotational_period_days = 27.322;
moon.radius_km = 1738;
moon.mu = 4.903e3;
moon.semi_major_axis_km = 3.844e5;
moon.orbital_period_days = 27.322;
moon.eccentriciy = .0549005;
moon.inclination_rad = 5.15*pi/180;

save -mat-binary moon.mat moon;

% mercury constants
mercury.axial_rotational_period_days = 58.6461;
mercury.radius_km = 2440.12;
mercury.mu = 2.203e4;
mercury.semi_major_axis_km = .387099*AU;
mercury.orbital_period_years = .2408;
mercury.eccentriciy = .205627;
mercury.inclination_rad = 7.00402*pi/180;

save -mat-binary mercury.mat mercury;

% venus constants
venus.axial_rotational_period_daysW = 243.01;
venus.radius_km = 6110;
venus.mu = 3.2486e5;
venus.semi_major_axis_km = .723332*AU;
venus.orbital_period_years = .6152;
venus.eccentriciy = .006793;
venus.inclination_rad = 3.39425*pi/180;

save -mat-binary venus.mat venus;

% earth constants
earth.axial_rotational_period_days = .99726;
earth.radius_km = 6378.14;
earth.mu = 3.986e5;
earth.semi_major_axis_km = AU;
earth.orbital_period_years = 1;
earth.eccentriciy = .016726;
earth.inclination_rad = 0;

save -mat-binary earth.mat earth;

% mars constants
mars.axial_rotational_period_days = 1.026;
mars.radius_km = 3394.74;
mars.mu = 4.2828e4;
mars.semi_major_axis_km = 1.523691*AU;
mars.orbital_period_years = 1.8808;
mars.eccentriciy = .093368;
mars.inclination_rad = 1.84992*pi/180;

save -mat-binary mars.mat mars;

% jupiter constants
jupiter.axial_rotational_period_days = .41354;
jupiter.radius_km = 70452;
jupiter.mu = 1.2671e8;
jupiter.semi_major_axis_km = 5.202803*AU;
jupiter.orbital_period_years = 11.86;
jupiter.eccentriciy = .048435;
jupiter.inclination_rad = 1.30618*pi/180;

save -mat-binary jupiter.mat jupiter;

% saturn constants
saturn.axial_rotational_period_days = .44403;
saturn.radius_km = 57822;
saturn.mu = 3.794e7;
saturn.semi_major_axis_km = 9.538843*AU;
saturn.orbital_period_years = 29.46;
saturn.eccentriciy = .055682;
saturn.inclination_rad = 2.48715*pi/180;

save -mat-binary saturn.mat saturn;

% uranus constants
uranus.axial_rotational_period_daysW = .68;
uranus.radius_km = 25150;
uranus.mu = 5.780e6;
uranus.semi_major_axis_km = 19.181951*AU;
uranus.orbital_period_years = 84;
uranus.eccentriciy = .047209;
uranus.inclination_rad = .77220*pi/180;

save -mat-binary uranus.mat uranus;

% neptune constants
neptune.axial_rotational_period_days = .57;
neptune.radius_km = 25092;
neptune.mu = 6.9e6;
neptune.semi_major_axis_km = 30.057779*AU;
neptune.orbital_period_years = 164.8;
neptune.eccentriciy = .008575;
neptune.inclination_rad = 1.77320*pi/180;

save -mat-binary neptune.mat neptune;

% pluto constants
pluto.axial_rotational_period_days = 6.3874;
pluto.radius_km = 2500;
pluto.mu = 1e3;
pluto.semi_major_axis_km = 39.4385*AU;
pluto.orbital_period_years = 247.7;
pluto.eccentriciy = .2481112;
pluto.inclination_rad = 17.16908*pi/180;

save -mat-binary pluto.mat pluto;

end