function [npoint, nline, nloop, nsurface] = set_Circles(fo, Ni, ...
    Diam, r, c, NACA, Pitch_D, skew, level, ref_angle, R, npoint, nline, nloop, nsurface, nres, ...
    n_airfoil, n_inlet, n_wake, n_thick)


mesh_param = 0.01;

%1 Adding points for the circles: first the center at each level, then 10
% points on the circle
fprintf(fo, '\n' );
fprintf(fo, 'Point(%d) = {%f, %f, %f, %d};\n', npoint + 1, 0, 0, r, mesh_param);
fprintf(fo, 'Point(%d) = {%f, %f, %f, %d};\n', npoint + 2, R * cos(ref_angle), R * sin(ref_angle), r, mesh_param);
fprintf(fo, 'Point(%d) = {%f, %f, %f, %d};\n', npoint + 3, R * cos(ref_angle + pi / 4), R * sin(ref_angle + pi / 4), r, mesh_param);
fprintf(fo, 'Point(%d) = {%f, %f, %f, %d};\n', npoint + 4, R * cos(ref_angle + pi / 2), R * sin(ref_angle + pi / 2), r, mesh_param);
fprintf(fo, 'Point(%d) = {%f, %f, %f, %d};\n', npoint + 5, R * cos(ref_angle + 6 * pi / 8), R * sin(ref_angle + 6 * pi / 8), r, mesh_param);
fprintf(fo, 'Point(%d) = {%f, %f, %f, %d};\n', npoint + 6, R * cos(ref_angle + 7 * pi / 8), R * sin(ref_angle + 7 * pi / 8), r, mesh_param);
fprintf(fo, 'Point(%d) = {%f, %f, %f, %d};\n', npoint + 7, R * cos(ref_angle + pi), R * sin(ref_angle + pi), r, mesh_param);
fprintf(fo, 'Point(%d) = {%f, %f, %f, %d};\n', npoint + 8, R * cos(ref_angle + 9 * pi / 8), R * sin(ref_angle + 9 * pi / 8), r, mesh_param);
fprintf(fo, 'Point(%d) = {%f, %f, %f, %d};\n', npoint + 9, R * cos(ref_angle + 10 * pi / 8), R * sin(ref_angle + 10 * pi / 8), r, mesh_param);
fprintf(fo, 'Point(%d) = {%f, %f, %f, %d};\n', npoint + 10, R * cos(ref_angle + 3 * pi / 2), R * sin(ref_angle + 3 * pi / 2), r, mesh_param);
fprintf(fo, 'Point(%d) = {%f, %f, %f, %d};\n', npoint + 11, R * cos(ref_angle - pi / 4), R * sin(ref_angle - pi / 4), r, mesh_param);


% points inside
R = 0.9 * R;
fprintf(fo, 'Point(%d) = {%f, %f, %f, %d};\n', npoint + 12, R * cos(ref_angle), R * sin(ref_angle), r, mesh_param);
fprintf(fo, 'Point(%d) = {%f, %f, %f, %d};\n', npoint + 13, R * cos(ref_angle + pi / 4), R * sin(ref_angle + pi / 4), r, mesh_param);
fprintf(fo, 'Point(%d) = {%f, %f, %f, %d};\n', npoint + 14, R * cos(ref_angle + pi / 2), R * sin(ref_angle + pi / 2), r, mesh_param);
fprintf(fo, 'Point(%d) = {%f, %f, %f, %d};\n', npoint + 15, R * cos(ref_angle + 6 * pi / 8), R * sin(ref_angle + 6 * pi / 8), r, mesh_param);
fprintf(fo, 'Point(%d) = {%f, %f, %f, %d};\n', npoint + 16, R * cos(ref_angle + 7 * pi / 8), R * sin(ref_angle + 7 * pi / 8), r, mesh_param);
fprintf(fo, 'Point(%d) = {%f, %f, %f, %d};\n', npoint + 17, R * cos(ref_angle + pi), R * sin(ref_angle + pi), r, mesh_param);
fprintf(fo, 'Point(%d) = {%f, %f, %f, %d};\n', npoint + 18, R * cos(ref_angle + 9 * pi / 8), R * sin(ref_angle + 9 * pi / 8), r, mesh_param);
fprintf(fo, 'Point(%d) = {%f, %f, %f, %d};\n', npoint + 19, R * cos(ref_angle + 10 * pi / 8), R * sin(ref_angle + 10 * pi / 8), r, mesh_param);
fprintf(fo, 'Point(%d) = {%f, %f, %f, %d};\n', npoint + 20, R * cos(ref_angle + 3 * pi / 2), R * sin(ref_angle + 3 * pi / 2), r, mesh_param);
fprintf(fo, 'Point(%d) = {%f, %f, %f, %d};\n', npoint + 21, R * cos(ref_angle - pi / 4), R * sin(ref_angle - pi / 4), r, mesh_param);



% add the circles
fprintf(fo, '\n' );
fprintf(fo, 'Circle(%d) = {%d, %d, %d};\n', nline + 1, npoint + 2, npoint + 1,  npoint + 3);
fprintf(fo, 'Circle(%d) = {%d, %d, %d};\n', nline + 2, npoint + 3, npoint + 1,  npoint + 4);
fprintf(fo, 'Circle(%d) = {%d, %d, %d};\n', nline + 3, npoint + 4, npoint + 1,  npoint + 5);
% fprintf(fo, 'Circle(%d) = {%d, %d, %d};\n', nline + 4, npoint + 5, npoint + 1,  npoint + 6);
% fprintf(fo, 'Circle(%d) = {%d, %d, %d};\n', nline + 5, npoint + 6, npoint + 1,  npoint + 7);
fprintf(fo, 'Circle(%d) = {%d, %d, %d};\n', nline + 4, npoint + 5, npoint + 1,  npoint + 7);
% fprintf(fo, 'Circle(%d) = {%d, %d, %d};\n', nline + 6, npoint + 7, npoint + 1,  npoint + 8);
% fprintf(fo, 'Circle(%d) = {%d, %d, %d};\n', nline + 7, npoint + 8, npoint + 1,  npoint + 9);
fprintf(fo, 'Circle(%d) = {%d, %d, %d};\n', nline + 6, npoint + 7, npoint + 1,  npoint + 9);
fprintf(fo, 'Circle(%d) = {%d, %d, %d};\n', nline + 8, npoint + 9, npoint + 1,  npoint + 10);
fprintf(fo, 'Circle(%d) = {%d, %d, %d};\n', nline + 9, npoint + 10, npoint + 1,  npoint + 11);
fprintf(fo, 'Circle(%d) = {%d, %d, %d};\n', nline + 10, npoint + 11, npoint + 1,  npoint + 2);




% Add some transfinite curves
fprintf(fo, '\n' );
fprintf(fo, 'Transfinite Curve {%d, %d} = %d Using Progression %f;\n', ...
    nline + 1, nline + 10, n_inlet, 1);
fprintf(fo, 'Transfinite Curve {%d, %d} = %d Using Progression %f;\n', ...
    nline + 2, nline + 9, n_airfoil, 1);
fprintf(fo, 'Transfinite Curve {%d, %d} = %d Using Progression %f;\n', ...
    nline + 3, nline + 8, n_wake, 1);
% fprintf(fo, 'Transfinite Curve {%d, %d} = %d Using Progression %f;\n', ...
%     nline + 4, nline + 7, n_extra, 1);
% fprintf(fo, 'Transfinite Curve {%d, %d} = %d Using Progression %f;\n', ...
%     nline + 5, nline + 3, n_thick, 1);
fprintf(fo, 'Transfinite Curve {%d, %d} = %d Using Progression %f;\n', ...
    nline + 4, nline + 6, n_thick, 1);




% 
%6 Update parameters

npoint = npoint + 21;
nline = nline + 10;
% nloop = nloop + 1;


end