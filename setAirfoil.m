

function [npoint, nline, nloop, nsurface] = setAirfoil(fo, Ni, ...
    Diam, r, r0, c, NACA, Pitch_D, skew, npoint, nline, nloop, nsurface, ...
    n_inlet, r_in_airfoil, n_airfoil)


% Extra parameters: the pitch and skew and the radius

% Diam = 10; % diameter
% r = 1; % current radius
% 
% Pitch_D = 0.6927; % 0.6927
% Pitch = Diam * Pitch_D; % Input pitch data
% alpha = atan(Pitch / (2 * pi * r)); % pitch angle
% theta_s = pi / 180 * 10.76; % 10.76 degree skew
% S = r * theta_s / cos(alpha); % Skew

Pitch = Diam * Pitch_D; % Input pitch data
alpha = atan(Pitch / (2 * pi * r0)); % pitch angle
theta_s = pi / 180 * skew; % 10.76 degree skew
S = r0 * theta_s / cos(alpha); % Skew



% Foil geometry
% c = 1;                      % Geometric chord length
s = 0.1;                    % Span (along y-axis)
aoa = 0;                    % Angle of attact (in degrees)
% inflation_factor_x = 5;     % Size of inflation layer w.r.t aifoil
inflation_factor_x = 20;     % Size of inflation layer w.r.t aifoil
% inflation_factor_z = 15;    % Size of inflation layer w.r.t aifoil
inflation_factor_z = 60;    % Size of inflation layer w.r.t aifoil

% NACA = [0 0 1 2];           % NACA 4-digit designation as a row vector;

% sizeNACA = 5e-3;  
sizeNACA = 1e-4;            % resolution at wall
sizeWake = 0.25; % 5e-2;            % resolution at wake
sizeFar  = 1;               % resolution at far field
shrink = 0.7;


% n_inlet = nres * 40;
% r_in_airfoil = 1.1;
% n_airfoil = 60;
xmax_c = 3.0;
xmin_c = -0.3;
ymax_c = 1.0*shrink;

x_center = 2;
y_center = 0;
R_circle = 4;



% xmax_c = 3.0;
% xmin_c = -0.3;
% ymax_c = 1.0*shrink;

% Surface resolution parameters
% Ni = 50;                  % Number of interpolation points along the foil

% ------------------------- END OF INPUT PARAMETER REGION -------------------- %


% ---------------------------------- LICENCE  -------------------------------- %
%                                                                              %
%     Copyrighted 2011, 2012 by Håkon Strandenes, hakostra@stud.ntnu.no        %
%                                                                              % 
%     This program is free software: you can redistribute it and/or modify     %
%     it under the terms of the GNU General Public License as published by     %
%     the Free Software Foundation, either version 3 of the License, or        %
%     (at your option) any later version.                                      %
%                                                                              %
%     This program is distributed in the hope that it will be useful,          %
%     but WITHOUT ANY WARRANTY; without even the implied warranty of           %
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            %
%     GNU General Public License for more details.                             %
%                                                                              %
%     You should have received a copy of the GNU General Public License        %
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.    %
% ---------------------------------------------------------------------------- %


% Create a vector with x-coordinates, camber and thickness
beta=linspace(0,pi,Ni);
x = c*(0.5*(1-cos(beta)));
z_c = zeros(size(x));
z_t = zeros(size(x));
theta = zeros(size(x));


% Values of m, p and t
m = NACA(1)/100;
p = NACA(2)/10;
t = (NACA(3)*10 + NACA(4))/100;


% Calculate thickness
% The upper expression will give the airfoil a finite thickness at the trailing
% edge, witch might cause trouble. The lower expression is corrected to give 
% zero thickness at the trailing edge, but the foil is strictly speaking no
% longer a proper NACA airfoil.
%
% See http://turbmodels.larc.nasa.gov/naca4412sep_val.html
%     http://en.wikipedia.org/wiki/NACA_airfoil

%z_t = (t*c/0.2) * (0.2969.*(x/c).^0.5 - 0.1260.*(x/c) - 0.3516.*(x/c).^2 + 0.2843.*(x/c).^3 - 0.1015.*(x/c).^4);
z_t = (t*c/0.2) * (0.2969.*(x/c).^0.5 - 0.1260.*(x/c) - 0.3516.*(x/c).^2 + 0.2843.*(x/c).^3 - 0.1036.*(x/c).^4);


% Calculate camber
if (p > 0)
  % Calculate camber
  z_c = z_c + (m.*x/p^2) .* (2*p - x/c) .* (x < p*c);
  z_c = z_c + (m.*(c-x)/(1-p)^2) .* (1 + x/c - 2*p) .* (x >= p*c);


  % Calculate theta-value
  theta = theta + atan( (m/p^2) * (2*p - 2*x/c) ) .* (x < p*c);
  theta = theta + atan( (m/(1-p)^2) * (-2*x/c + 2*p) ) .* (x >= p*c);
end


% Calculate coordinates of upper surface
Xu = x - z_t.*sin(theta);
Zu = z_c + z_t.*cos(theta);

% Calculate coordinates of lower surface
Xl = x + z_t.*sin(theta);
Zl = z_c - z_t.*cos(theta);

upper = [Xu ; Zu];
lower = [Xl ; Zl];


% Merge upper and lower surface (NB: Assume that the trailing edge is sharp)
% (see comments w.r.t. thickness calculation above)
X = [ upper(1,:) lower(1,Ni-1:-1:2) ];
Z = [ upper(2,:) lower(2,Ni-1:-1:2) ];

Rot = [cos(alpha), sin(alpha); - sin(alpha), cos(alpha)];

% after or before adding the skew ?
V = [0.5 * c + S; 0] + Rot * [X - 0.5 * c; Z]; 
X = V(1, :);
Z = V(2, :);

Vect = [0.5 * c + S; 0] + Rot * ([xmax_c, xmax_c, xmin_c, xmin_c, 1, 1, xmax_c, -sqrt(ymax_c^2+xmin_c^2); ...
    ymax_c * 1.25, -ymax_c * 1.25, ymax_c, -ymax_c, ymax_c, -ymax_c, 0, 0] - [0.5 * c; 0]); 

N = length(X);

% Open file
% fo = fopen(fileName, 'w');

fprintf(fo, 'aoa = %f;\n', aoa);
fprintf(fo, 'sizeNACA = %f;\n', sizeNACA);
fprintf(fo, 'sizeWake = %f;\n', sizeWake);
fprintf(fo, 'sizeFar = %f;\n', sizeFar);

% Add some constants to print

fprintf(fo, 'shrink = %f;\n', shrink);
fprintf(fo, 'xmax_c = %f;\n', 3.0);
fprintf(fo, 'xmin_c = %f;\n', -0.3);
fprintf(fo, 'ymax_c = %f;\n', 1.0*shrink);
fprintf(fo, 'n_inlet = %d;\n', 40);
fprintf(fo, 'r_inlet = %f;\n', 0.3);
fprintf(fo, 'r_in_airfoil = %f;\n', 1.1);
fprintf(fo, 'n_inlet = %d;\n', 120*shrink+50);
fprintf(fo, 'r_vertical = %f;\n', 1.015);
fprintf(fo, 'n_airfoil = %d;\n', 170);
fprintf(fo, 'n_wake = %d;\n', 220);
fprintf(fo, 'r_wake = %f;\n', 1.001);
fprintf(fo, 'r_wake_airfoil = %f;\n', 1.005);
fprintf(fo, 'r_vertical_wake = %f;\n', 1.0);
fprintf(fo, 'x_center = %f;\n', x_center);
fprintf(fo, 'y_center = %f;\n', y_center);
fprintf(fo, 'R_circle = %f;\n', R_circle);


fprintf(fo, '\n' );

%1 Write points
for i=1:N
	disp(i)
	fprintf(fo, 'Point(%d) = {%f, %f, %f, sizeNACA};\n', npoint + i, X(i), Z(i), r);
end

fprintf(fo, '\n' );

fprintf(fo, 'Spline(%d) = {', nline + 1);
for i = 1:Ni-1
    fprintf(fo, '%d, ', npoint + i);
end
fprintf(fo, '%d', npoint + Ni);
fprintf(fo, '};\n' );

fprintf(fo, 'Spline(%d) = {', nline + 2);
for i = 1:Ni-1
    fprintf(fo, '%d, ', npoint + Ni+i-1);
end
fprintf(fo, '%d',  npoint + 1);
fprintf(fo, '};\n' );


fprintf(fo, '\n' );
Nsplit = round(0.24 * Ni);
fprintf(fo, 'Split Curve {%d} Point {%d};\n', nline + 1, npoint + Nsplit);
fprintf(fo, 'Split Curve {%d} Point {%d};\n', nline + 2, npoint + 2 * Ni - Nsplit);


%3 Add some transfinite curves

fprintf(fo, '\n' );
fprintf(fo, 'Transfinite Curve {%d, %d} = %d Using Bump %f;\n', ...
    nline + 6, - (nline + 3), n_inlet, r_in_airfoil);
fprintf(fo, 'Transfinite Curve {%d, %d} = %d Using Progression %f;\n', ...
    nline + 5, nline + 4, n_airfoil, 1);

%4 Update parameters
npoint = npoint + N;
nline = nline + 6;


end

