clear all
close all


fileName='blade.geo';

% Open file
fo = fopen(fileName, 'w');


npoint = 0;
nline = 0;
nloop = 0;
nsurface = 0;

%% Blade data

Diam = 0.08;
Rmax = Diam / 2;

radiusq = [0.16; 0.25; 0.30; 0.40; 0.50; 0.60; 0.70; 0.80; 0.90; 0.95] * Rmax;
Pitches_Dq = [0.5783; 0.6130; 0.6310; 0.66631; 0.6916; 0.7121; 0.7213; 0.7161; 0.6928; 0.6749];
Cq = [0.1529; 0.1772; 0.1892; 0.2093; 0.2247; 0.2335; 0.2338; 0.2191; 0.1808; 0.1420] * Diam;

t0_D = [46.61; 42.15; 38.50; 32.00; 26.00; 20.54; 15.59; 11.04; 6.99; 4.74];
NACA3s = floor(100 * t0_D * Diam ./ Cq / 1000 / 10);
NACA4s = 100 * t0_D * Diam ./ Cq / 1000 - 10 * NACA3s;

NACAsq = [[3.16 0 1 2]; ...
         [3.49 0 1 2]; ...
         [3.57 0 1 2]; ...
         [3.38 0 1 2]; ...
         [2.94 0 1 2]; ...
         [2.50 0 1 2]; ...
         [2.19 0 1 2]; ...
         [1.98 0 1 2]; ...
         [1.62 0 1 2]; ...
         [1.28 0 1 2]];
NACAsq(:, 3) = NACA3s;
NACAsq(:, 4) = NACA4s;
Skewsq = [-2.63; -3.99; -4.39; -4.39; -3.14; -0.82; 2.49; 6.33; 10.64; 12.93];



%% Discretization parameters

res = 0.25;
Nz = res * 60;
Ni = 70; % Number of points used to draw half an airfoil
nvert = ceil(100 * res);% number of discretizaton points between two levels vertically
n_thick = ceil(100 * res); % number of discretizaton points along the thickness of the wake in the Cmesh
n_inlet = ceil(res * 40); % number of discretizaton points along the top circular parts of the blade
n_wake = ceil(res * 170);% number of discretizaton points along the wake in the Cmesh
r_in_airfoil = 1.1; % transfinite bump progression along the top circular parts of the blade
n_airfoil_lateral = ceil(10 * res); % number of discretization points along the line that cuts the bottom and top airfoil in half
n_airfoil = ceil(130 * res);% number of discretizaton points along the bottom long parts of the airfoil
n_splines = ceil(150 * res); % number of discretization points along the splines
r_splines = 1.05; % transfinite progression parameter along the splines
R = 10 * max(Cq); % the radius of the cylinder


%% Domain extension 
Pitches_Dq = [Pitches_Dq(1); Pitches_Dq; Pitches_Dq(end)];
Cq = [(2 * Cq(1) - Cq(2)); Cq; (2 * Cq(end) - Cq(end - 1))];
t0_D = [t0_D(1); t0_D; t0_D(end)];
NACAsq = [NACAsq(1, :); NACAsq;  NACAsq(end, :)];
Skewsq = [Skewsq(1); Skewsq; Skewsq(end)];
N0 = size(radiusq, 1);
rmin = radiusq(1);
rmax = radiusq(end);
radiusq = [-2 * Rmax; radiusq; 3 * Rmax];

radius0 = linspace(radiusq(1), rmin, Nz).';
radius1 = linspace(rmin, rmax, Nz).';
radius2 = linspace(rmax, radiusq(end), Nz).';
radius = [radius0(1 : end - 1); radius1; radius2(2 : end)];
level_min = length(radius0);
level_max = level_min + length(radius1) - 1;

Pitches_D = pchip(radiusq,Pitches_Dq,radius);
C = pchip(radiusq,Cq,radius);
for i = 1 : 4
    NACAs(:, i) = pchip(radiusq,NACAsq(:, i),radius);
end
Skews = pchip(radiusq,Skewsq,radius);




%% Setting the airfoils

nradii = length(radius);
% nradii = 2;
% Loop over radii
for i = 1 : nradii
% for i = 1 : 2
   r = radius(i);
   if(r < rmin)
       r0 = rmin;
%        nvert = nvert2;
   elseif(r > rmax)
       r0 = rmax;
%        nvert = nvert2;
   else
       r0 = r;
%        nvert = nvert1;
   end   
   Pitch_D = Pitches_D(i);
   skew = Skews(i);
   c = C(i);
   NACA = NACAs(i, :);
   [npoint, nline, nloop, nsurface] = setAirfoil(fo, Ni, Diam, r, r0, c, ...
    NACA, Pitch_D, skew, npoint, nline, nloop, nsurface, n_inlet, ...
    r_in_airfoil, n_airfoil);
   if(i == 1)
       npoint_per_level = npoint;
   end
end
npoint_airfoil = npoint;
nline_airfoil = nline;



%% The contours of the Cmeshes

for i = 1 : nradii
% for i = 1 : 2
   r = radius(i);
   if(r < rmin)
       r0 = rmin;
%        nvert = nvert2;
   elseif(r > rmax)
       r0 = rmax;
%        nvert = nvert2;
   else
       r0 = r;
%        nvert = nvert1;
   end  
   Pitch_D = Pitches_D(i);
   skew = Skews(i);
   c = C(i);
   NACA = NACAs(i, :);
    [npoint, nline, nloop, nsurface] = setC_Circle_contours(fo, Ni, ...
        Diam, r, r0, c, NACA, Pitch_D, skew, i, npoint, ...
        nline, nloop, nsurface, res, n_airfoil, n_inlet, n_wake, n_thick);
end

npoint_Cmesh = npoint;
nline_Cmesh = nline;



%% Add the cylinder around the C-mesh contours

Pitches = Diam * Pitches_Dq; % Input pitch data
alphas = atan(Pitches ./ (2 * pi * radiusq)); % pitch angle
% ref_angle = pi + mean(alphas)

% ref_angle = 3 * pi / 4;
ref_angle = pi;

% for i = 1 : 2
for i = 1 : nradii
   r = radius(i);
   Pitch_D = Pitches_D(i);
   skew = Skews(i);
   c = C(i);
   NACA = NACAs(i, :);
   if(r < rmin)
       r0 = rmin;
   elseif(r > rmax)
       r0 = rmax;
   else
       r0 = r;
   end     
    [npoint, nline, nloop, nsurface] = set_Circles(fo, Ni, Diam, r, c, ...
        NACA, Pitch_D, skew, i, ref_angle, R, npoint, nline, nloop, ...
        nsurface, res, n_airfoil, n_inlet, n_wake, n_thick);
end



%% Set the Cmeshes and inter-volumes


[npoint, nline, nloop, nsurface] = setCmeshes(fo, Ni, npoint, nline, ...
    nloop, nsurface, npoint_per_level, nline_airfoil, npoint_airfoil, ...
    nradii, res, n_airfoil, n_inlet, n_wake, nvert, n_thick);


%% Set the exterior meshes and intervolumes

[npoint, nline, nloop, nsurface] = setExtmeshes(fo, Ni, npoint, nline, ...
    nloop, nsurface, npoint_per_level, nline_airfoil, npoint_airfoil, ...
    npoint_Cmesh, nline_Cmesh, nradii, res, n_airfoil, n_inlet, n_wake, ...
    nvert, n_thick, n_splines, r_splines);



%% Set the extensions top and bottom

[npoint, nline, nloop, nsurface] = setBladeVolumesExtension(fo, Ni, ...
    npoint, nline, nloop, nsurface, nradii, res, n_airfoil_lateral, ...
    radius, rmin, rmax);


%% Set the physical boundaries

setPhysicalBoundaries(fo, radius, level_min, level_max);



fclose(fo);