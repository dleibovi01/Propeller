function [npoint, nline, nloop, nsurface] = setC_Circle_contours(fo, Ni, ...
    Diam, r, r0, c, NACA, Pitch_D, skew, level, npoint, nline, nloop, nsurface, nres, ...
    n_airfoil, n_inlet, n_wake, n_thick)



mesh_param = 0.01;


Pitch = Diam * Pitch_D; % Input pitch data
alpha = atan(Pitch / (2 * pi * r0)); % pitch angle
theta_s = pi / 180 * skew; % 10.76 degree skew
S = r0 * theta_s / cos(alpha); % Skew

% Foil geometry
% c = 1;                      % Geometric chord length
s = 0.1;                    % Span (along y-axis)
aoa = 0;                    % Angle of attact (in degrees)
inflation_factor_x = 20;     % Size of inflation layer w.r.t aifoil
inflation_factor_z = 60;    % Size of inflation layer w.r.t aifoil

sizeNACA = 5e-3;            % resolution at wall
sizeWake = 0.25; % 5e-2;            % resolution at wake
sizeFar  = 1;               % resolution at far field
shrink = 0.7;




% n_inlet = nres * 40;
r_inlet = 0.3;
n_vertical = nres * (60*shrink+25);
% n_airfoil = 60;
% n_wake = nres * 110;
r_wake = 1.001;
r_vertical_wake = 1.0;
xmax_c = 3.0 * c;
xmin_c = -0.3 * c;
ymax_c = 1.0*shrink * c;

% x_center = 1 * c + S;
% y_center = - 0.5 * c;
x_center = 0.01;
y_center = 0;



% x_center = S +  c * cos(alpha);
% y_center = - c * sin(alpha);

% n_circle = nres * 120;



% Values of m, p and t
m = NACA(1)/100;
p = NACA(2)/10;
t = (NACA(3)*10 + NACA(4))/100;


Rot = [cos(alpha), sin(alpha); - sin(alpha), cos(alpha)];
Vect = [0.5 * c + S; 0] + Rot * ([xmax_c, xmax_c, xmin_c, xmin_c, c, c, xmax_c, -sqrt(ymax_c^2+xmin_c^2); ...
    ymax_c * 1.25, -ymax_c * 1.25, ymax_c, -ymax_c, ymax_c, -ymax_c, 0, 0] - [0.5 * c; 0]); 

Vect_ext = [0.5 * c + S; 0] + Rot * 1.1 * ([xmax_c, xmax_c, xmin_c, xmin_c, c, c, xmax_c, -sqrt(ymax_c^2+xmin_c^2); ...
    ymax_c * 1.25, -ymax_c * 1.25, ymax_c, -ymax_c, ymax_c, -ymax_c, 0, 0] - [0.5 * c; 0]); 



%1 Add some points for the exterior of the C-mesh
fprintf(fo, '\n' );
fprintf(fo, 'Point(%d) = {%f, %f, %f, %d};\n', npoint + 1, Vect(1, 1), Vect(2, 1), r, mesh_param);
fprintf(fo, 'Point(%d) = {%f, %f, %f, %d};\n', npoint + 2, Vect(1, 2), Vect(2, 2), r, mesh_param);
fprintf(fo, 'Point(%d) = {%f, %f, %f, %d};\n', npoint + 3, Vect(1, 3), Vect(2, 3), r, mesh_param);
fprintf(fo, 'Point(%d) = {%f, %f, %f, %d};\n', npoint + 4, Vect(1, 4), Vect(2, 4), r, mesh_param);
fprintf(fo, 'Point(%d) = {%f, %f, %f, %d};\n', npoint + 5, Vect(1, 5), Vect(2, 5), r, mesh_param);
fprintf(fo, 'Point(%d) = {%f, %f, %f, %d};\n', npoint + 6, Vect(1, 6), Vect(2, 6), r, mesh_param);
fprintf(fo, 'Point(%d) = {%f, %f, %f, %d};\n', npoint + 7, Vect(1, 7), Vect(2, 7), r, mesh_param);
fprintf(fo, 'Point(%d) = {%f, %f, %f, %d};\n', npoint + 8, Vect(1, 8), Vect(2, 8), r, mesh_param);

% Add points exterior to the C-mesh

Nsplit = round(0.24 * Ni);
fprintf(fo, '\n' );

fprintf(fo, 'xyz_%d[] = Point{%d};\n', (level - 1) * 2 * (Ni - 1) + 1, (level - 1) * 2 * (Ni - 1) + 1);
fprintf(fo, 'Point(%d) = {%f + 1 * (%f - xyz_%d[0]), %f + 1 * (%f - xyz_%d[1]), %f, %d};\n', ...
    npoint + 9, Vect(1, 8), Vect(1, 8), (level - 1) * 2 * (Ni - 1) + 1, Vect(2, 8), Vect(2, 8), (level - 1) * 2 * (Ni - 1) + 1, r, mesh_param);

fprintf(fo, 'xyz_%d[] = Point{%d};\n', (level - 1) * 2 * (Ni - 1) + 2, (level - 1) * 2 * (Ni - 1) + Nsplit);
fprintf(fo, 'Point(%d) = {%f + 1 * (%f - xyz_%d[0]), %f + 1 * (%f - xyz_%d[1]), %f, %d};\n', ...
    npoint + 10, Vect(1, 3), Vect(1, 3), (level - 1) * 2 * (Ni - 1) + 2, Vect(2, 3), Vect(2, 3), (level - 1) * 2 * (Ni - 1) + 2, r, mesh_param);

fprintf(fo, 'xyz_%d[] = Point{%d};\n', (level - 1) * 2 * (Ni - 1) + 3, (level - 1) * 2 * (Ni - 1) + Ni);
fprintf(fo, 'Point(%d) = {%f + 1 * (%f - xyz_%d[0]), %f + 1 * (%f - xyz_%d[1]), %f, %d};\n', ...
    npoint + 11, Vect(1, 5), Vect(1, 5), (level - 1) * 2 * (Ni - 1) + 3, Vect(2, 5), Vect(2, 5), (level - 1) * 2 * (Ni - 1) + 3, r, mesh_param);

% This one to be modified
% fprintf(fo, 'Point(%d) = {%f + 1 * (%f - %f), %f + 1 * (%f - %f), %f, %d};\n', ...
%     npoint + 12, Vect(1, 1), Vect(1, 1), Vect(1, 7), Vect(2, 1), Vect(2, 1), Vect(2, 7), r, mesh_param);
% fprintf(fo, 'Point(%d) = {2 * %f - 0.5 * %f - 0.25 * %f, 2 * %f - 0.5 * %f - 0.25 * %f, %f, %d};\n', ...
%     npoint + 12, Vect(1, 1), Vect(1, 7), Vect(1, 5), Vect(2, 1), Vect(2, 7), Vect(2, 5), r, mesh_param);
% fprintf(fo, 'Point(%d) = {3 * %f - 1 * %f - 0.5 * %f, 3 * %f - 1 * %f - 0.5 * %f, %f, %d};\n', ...
%     npoint + 12, Vect(1, 1), Vect(1, 7), Vect(1, 5), Vect(2, 1), Vect(2, 7), Vect(2, 5), r, mesh_param);
fprintf(fo, 'Point(%d) = {4 * %f - 2 * %f - 1 * %f, 4 * %f - 2 * %f - 1 * %f, %f, %d};\n', ...
    npoint + 12, Vect(1, 1), Vect(1, 7), Vect(1, 5), Vect(2, 1), Vect(2, 7), Vect(2, 5), r, mesh_param);


fprintf(fo, 'Point(%d) = {%f + 0.5 * (%f - %f), %f + 0.5 * (%f - %f), %f, %d};\n', ...
    npoint + 13, Vect(1, 1), Vect(1, 1), Vect(1, 5), Vect(2, 1), Vect(2, 1), Vect(2, 5), r, mesh_param);


fprintf(fo, 'xyz_%d[] = Point{%d};\n', (level - 1) * 2 * (Ni - 1) + 4, (level - 1) * 2 * (Ni - 1) + Ni);
fprintf(fo, 'Point(%d) = {%f + 1 * (%f - xyz_%d[0]), %f + 1 * (%f - xyz_%d[1]), %f, %d};\n', ...
    npoint + 14, Vect(1, 7), Vect(1, 7), (level - 1) * 2 * (Ni - 1) + 4, Vect(2, 7), Vect(2, 7), (level - 1) * 2 * (Ni - 1) + 4, r, mesh_param);

fprintf(fo, 'Point(%d) = {%f + 0.5 * (%f - %f), %f + 0.5 * (%f - %f), %f, %d};\n', ...
    npoint + 15, Vect(1, 2), Vect(1, 2), Vect(1, 6), Vect(2, 2), Vect(2, 2), Vect(2, 6), r, mesh_param);

% This one to be modified
% fprintf(fo, 'Point(%d) = {%f + 1 * (%f - %f), %f + 1 * (%f - %f), %f, %d};\n', ...
%     npoint + 16, Vect(1, 2), Vect(1, 2), Vect(1, 7), Vect(2, 2), Vect(2, 2), Vect(2, 7), r, mesh_param);
fprintf(fo, 'Point(%d) = {4 * %f - 2 * %f - 1 * %f, 4 * %f - 2 * %f - 1 * %f, %f, %d};\n', ...
    npoint + 16, Vect(1, 2), Vect(1, 7), Vect(1, 6), Vect(2, 2), Vect(2, 7), Vect(2, 6), r, mesh_param);

fprintf(fo, 'xyz_%d[] = Point{%d};\n', (level - 1) * 2 * (Ni - 1) + 5, (level - 1) * 2 * (Ni - 1) + Ni);
fprintf(fo, 'Point(%d) = {%f + 1 * (%f - xyz_%d[0]), %f + 1 * (%f - xyz_%d[1]), %f, %d};\n', ...
    npoint + 17, Vect(1, 6), Vect(1, 6), (level - 1) * 2 * (Ni - 1) + 3, Vect(2, 6), Vect(2, 6), (level - 1) * 2 * (Ni - 1) + 5, r, mesh_param);

fprintf(fo, 'xyz_%d[] = Point{%d};\n', (level - 1) * 2 * (Ni - 1) + 6, (level - 1) * 2 * (Ni - 1) + 2 * Ni - Nsplit);
fprintf(fo, 'Point(%d) = {%f + 1 * (%f - xyz_%d[0]), %f + 1 * (%f - xyz_%d[1]), %f, %d};\n', ...
    npoint + 18, Vect(1, 4), Vect(1, 4), (level - 1) * 2 * (Ni - 1) + 3, Vect(2, 4), Vect(2, 4), (level - 1) * 2 * (Ni - 1) + 6, r, mesh_param);



%2 Add some lines
fprintf(fo, '\n' );
Nsplit = round(0.24 * Ni);

fprintf(fo, 'Line(%d) = {%d, %d};\n', nline + 1, npoint + 3, npoint + 5);
fprintf(fo, 'Line(%d) = {%d, %d};\n', nline + 2, npoint + 4, npoint + 6);
fprintf(fo, 'Line(%d) = {%d, %d};\n', nline + 3, npoint + 5, npoint + 1);
fprintf(fo, 'Line(%d) = {%d, %d};\n', nline + 4, npoint + 6, npoint + 2);
fprintf(fo, 'Line(%d) = {%d, %d};\n', nline + 5, npoint + 1, npoint + 7);
fprintf(fo, 'Line(%d) = {%d, %d};\n', nline + 6, npoint + 2, npoint + 7);



%3 Add some circles
fprintf(fo, '\n' );
fprintf(fo, 'Circle(%d) = {%d, %d, %d};\n', nline + 7, npoint + 4,  (level - 1) * 2 * (Ni - 1) + 1, npoint + 8);
fprintf(fo, 'Circle(%d) = {%d, %d, %d};\n', nline + 8, npoint + 3, (level - 1) * 2 * (Ni - 1) + 1, npoint + 8);


%4 Add some transfinite curves
fprintf(fo, '\n' );
fprintf(fo, 'Transfinite Curve {%d, %d} = %d Using Progression %f;\n', ...
    nline + 3, nline + 4, n_wake, 1);
fprintf(fo, 'Transfinite Curve {%d, %d} = %d Using Progression %f;\n', ...
    -(nline + 5), -(nline + 6), n_thick, r_wake);
% fprintf(fo, 'Transfinite Curve {%d, %d} = %d Using Progression %f;\n', ...
%     - (nline + 7), - (nline + 8), n_inlet, r_inlet);
fprintf(fo, 'Transfinite Curve {%d, %d} = %d Using Progression %f;\n', ...
    (nline + 7), (nline + 8), n_inlet, 1);
fprintf(fo, 'Transfinite Curve {%d, %d} = %d Using Progression %f;\n', ...
    nline + 1, nline + 2, n_airfoil, 1);


% fprintf(fo, 'Transfinite Curve {%d, %d, %d, %d} = %d Using Progression %f;\n', ...
%     nline + 15, nline + 16, nline + 17, nline + 18, n_circle, 1);



%5 Add some Curved loops

fprintf(fo, '\n' );
fprintf(fo, 'Curve Loop(%d) = {%d, %d, %d, %d, %d, %d, %d, %d};\n', ...
    nloop + 1, nline + 3, nline + 5, nline + 7, - (nline + 8), ... 
    - (nline + 6), -(nline + 4), nline + 13, -(nline + 14));
% fprintf(fo, 'Curve Loop(%d) = {%d, %d, %d, %d};\n', ...
%     nloop + 2, nline + 15, nline + 16, nline + 17, nline + 18);

%5 Add some surfaces

% fprintf(fo, '\n' );
% fprintf(fo, 'Plane Surface(%d) = {%d, %d};\n', nsurface + 1, nloop + 1, ...
%     nloop + 2);
% 
% fprintf(fo, '\n' );
% fprintf(fo, 'Recombine Surface {%d};\n', nsurface + 1);
% 
%6 Update parameters

npoint = npoint + 18;
nline = nline + 8;
nloop = nloop + 1;


end