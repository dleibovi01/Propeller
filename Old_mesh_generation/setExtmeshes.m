

function [npoint, nline, nloop, nsurface] = setExtmeshes(fo, Ni, ...
    npoint, nline, nloop, nsurface, ...
    npoint_per_level, nline_airfoil, npoint_airfoil, npoint_Cmesh, ...
    nline_Cmesh, nradii, nres, n_airfoil, n_inlet, n_wake, n_vert, ...
    n_thick, n_splines, r_splines)



% NACA = [0 0 1 2];           % NACA 4-digit designation as a row vector;

sizeNACA = 5e-3;            % resolution at wall
sizeWake = 0.25; % 5e-2;            % resolution at wake
sizeFar  = 1;               % resolution at far field
shrink = 0.7;



%% Set the splines between the Cmesh and the circles, and the vertical lines
% between the circles

for i = 1 : nradii
% for i = 1 : 2
    fprintf(fo, '\n' );
    % Create the horizontal lines 
    nline_per_level = 20;
    npoint_C_level = 18;
    npoint_circle_level = 21;


    fprintf(fo, 'L_ext_%d = newc;\n', nline_per_level * (i - 1) + 1);
    fprintf(fo, 'BSpline(L_ext_%d) = {%d, %d, %d, %d};\n', ...
        nline_per_level * (i - 1) + 1, ...
        npoint_airfoil + (i - 1) * npoint_C_level + 3, ...
        npoint_airfoil + (i - 1) * npoint_C_level + 10, ...
        npoint_Cmesh + (i - 1) * npoint_circle_level + 21, ...
        npoint_Cmesh + (i - 1) * npoint_circle_level + 11);

    fprintf(fo, 'L_ext_%d = newc;\n', nline_per_level * (i - 1) + 2);
    fprintf(fo, 'BSpline(L_ext_%d) = {%d, %d, %d, %d};\n', ...
        nline_per_level * (i - 1) + 2, ...
        npoint_airfoil + (i - 1) * npoint_C_level + 5, ...
        npoint_airfoil + (i - 1) * npoint_C_level + 11, ...
        npoint_Cmesh + (i - 1) * npoint_circle_level + 20, ...
        npoint_Cmesh + (i - 1) * npoint_circle_level + 10);    
    
    fprintf(fo, 'L_ext_%d = newc;\n', nline_per_level * (i - 1) + 3);
    fprintf(fo, 'BSpline(L_ext_%d) = {%d, %d, %d, %d};\n', ...
        nline_per_level * (i - 1) + 3, ...
        npoint_airfoil + (i - 1) * npoint_C_level + 1, ...
        npoint_airfoil + (i - 1) * npoint_C_level + 12, ...
        npoint_Cmesh + (i - 1) * npoint_circle_level + 19, ...
        npoint_Cmesh + (i - 1) * npoint_circle_level + 9);      
  
%     fprintf(fo, 'L_ext_%d = newc;\n', nline_per_level * (i - 1) + 4);
%     fprintf(fo, 'BSpline(L_ext_%d) = {%d, %d, %d, %d};\n', ...
%         nline_per_level * (i - 1) + 4, ...
%         npoint_airfoil + (i - 1) * npoint_C_level + 1, ...
%         npoint_airfoil + (i - 1) * npoint_C_level + 13, ...
%         npoint_Cmesh + (i - 1) * npoint_circle_level + 18, ...
%         npoint_Cmesh + (i - 1) * npoint_circle_level + 8);     
     
    fprintf(fo, 'L_ext_%d = newc;\n', nline_per_level * (i - 1) + 5);
    fprintf(fo, 'BSpline(L_ext_%d) = {%d, %d, %d, %d};\n', ...
        nline_per_level * (i - 1) + 5, ...
        npoint_airfoil + (i - 1) * npoint_C_level + 7, ...
        npoint_airfoil + (i - 1) * npoint_C_level + 14, ...
        npoint_Cmesh + (i - 1) * npoint_circle_level + 17, ...
        npoint_Cmesh + (i - 1) * npoint_circle_level + 7);       
    
%     fprintf(fo, 'L_ext_%d = newc;\n', nline_per_level * (i - 1) + 6);
%     fprintf(fo, 'BSpline(L_ext_%d) = {%d, %d, %d, %d};\n', ...
%         nline_per_level * (i - 1) + 6, ...
%         npoint_airfoil + (i - 1) * npoint_C_level + 2, ...
%         npoint_airfoil + (i - 1) * npoint_C_level + 15, ...
%         npoint_Cmesh + (i - 1) * npoint_circle_level + 16, ...
%         npoint_Cmesh + (i - 1) * npoint_circle_level + 6);    
    
    fprintf(fo, 'L_ext_%d = newc;\n', nline_per_level * (i - 1) + 7);
    fprintf(fo, 'BSpline(L_ext_%d) = {%d, %d, %d, %d};\n', ...
        nline_per_level * (i - 1) + 7, ...
        npoint_airfoil + (i - 1) * npoint_C_level + 2, ...
        npoint_airfoil + (i - 1) * npoint_C_level + 16, ...
        npoint_Cmesh + (i - 1) * npoint_circle_level + 15, ...
        npoint_Cmesh + (i - 1) * npoint_circle_level + 5);       
    
    fprintf(fo, 'L_ext_%d = newc;\n', nline_per_level * (i - 1) + 8);
    fprintf(fo, 'BSpline(L_ext_%d) = {%d, %d, %d, %d};\n', ...
        nline_per_level * (i - 1) + 8, ...
        npoint_airfoil + (i - 1) * npoint_C_level + 6, ...
        npoint_airfoil + (i - 1) * npoint_C_level + 17, ...
        npoint_Cmesh + (i - 1) * npoint_circle_level + 14, ...
        npoint_Cmesh + (i - 1) * npoint_circle_level + 4);        
    
    fprintf(fo, 'L_ext_%d = newc;\n', nline_per_level * (i - 1) + 9);
    fprintf(fo, 'BSpline(L_ext_%d) = {%d, %d, %d, %d};\n', ...
        nline_per_level * (i - 1) + 9, ...
        npoint_airfoil + (i - 1) * npoint_C_level + 4, ...
        npoint_airfoil + (i - 1) * npoint_C_level + 18, ...
        npoint_Cmesh + (i - 1) * npoint_circle_level + 13, ...
        npoint_Cmesh + (i - 1) * npoint_circle_level + 3);    
    
    fprintf(fo, 'L_ext_%d = newc;\n', nline_per_level * (i - 1) + 10);
    fprintf(fo, 'BSpline(L_ext_%d) = {%d, %d, %d, %d};\n', ...
        nline_per_level * (i - 1) + 10, ...
        npoint_airfoil + (i - 1) * npoint_C_level + 8, ...
        npoint_airfoil + (i - 1) * npoint_C_level + 9, ...
        npoint_Cmesh + (i - 1) * npoint_circle_level + 12, ...
        npoint_Cmesh + (i - 1) * npoint_circle_level + 2);  
    
    % And now the vertical lines
    if(i < nradii)
        for j = 0 : 9
            if( j~= 3 && j ~= 5)
                fprintf(fo, 'L_ext_%d = newc;\n', nline_per_level * (i - 1) + 11 + j);
                fprintf(fo, 'Line(L_ext_%d) = {%d, %d};\n', ...
                    nline_per_level * (i - 1) + 11 + j, ...
                    npoint_Cmesh + (i - 1) * npoint_circle_level + 11 - j, ...
                    npoint_Cmesh + i * npoint_circle_level + 11 - j);
            end
        end
        
    end
    
end


%% Make the splines and the vertical lines transfinite

for i = 1 : nradii
    fprintf(fo, '\n' );
    for j = 1 : 10
        if(j ~= 4 && j~= 6)
            fprintf(fo, 'Transfinite Curve {L_ext_%d} = %d Using Progression %f;\n', ...
                nline_per_level * (i - 1) + j, n_splines, r_splines);
        end
    end  
    if(i < nradii)
        for j = 1 : 10
            if(j ~= 4 && j~= 6)
                fprintf(fo, 'Transfinite Curve {L_ext_%d} = %d Using Progression %f;\n', ...
                    nline_per_level * (i - 1) + 10 + j, n_vert, 1);    
            end
        end
    end
end


%% Declare the curve loops

for i = 1 : nradii - 1
    % newcl
    nloop_per_level = 30;
    n_ext_line_per_level = 8;
    nline_airfoil_per_level = 18;
    ncirc_line_per_level = 10;
    fprintf(fo, '\n' );
    
    
    % Horizontal faces    
    
    fprintf(fo, 'CL_ext_%d = newcl;\n', nloop_per_level * (i - 1) + 1);
    fprintf(fo, 'Curve Loop(CL_ext_%d) = {- %d, L_ext_%d, - %d, - L_ext_%d};\n', ...
        nloop_per_level * (i - 1) + 1, ...
        nline_airfoil + (i - 1) * n_ext_line_per_level + 1, ...
        nline_per_level * (i - 1) + 1, ...
        nline_Cmesh + (i - 1) * ncirc_line_per_level + 9, ...
        nline_per_level * (i - 1) + 2);      
    
    fprintf(fo, 'CL_ext_%d = newcl;\n', nloop_per_level * (i - 1) + 2);
    fprintf(fo, 'Curve Loop(CL_ext_%d) = {- %d, L_ext_%d, - %d, - L_ext_%d};\n', ...
        nloop_per_level * (i - 1) + 2, ...
        nline_airfoil + (i - 1) * n_ext_line_per_level + 3, ...
        nline_per_level * (i - 1) + 2, ...
        nline_Cmesh + (i - 1) * ncirc_line_per_level + 8, ...
        nline_per_level * (i - 1) + 3);    
    
    fprintf(fo, 'CL_ext_%d = newcl;\n', nloop_per_level * (i - 1) + 3);
    fprintf(fo, 'Curve Loop(CL_ext_%d) = {- %d, L_ext_%d, - %d, - L_ext_%d};\n', ...
        nloop_per_level * (i - 1) + 3, ...
        nline_airfoil + (i - 1) * n_ext_line_per_level + 5, ...
        nline_per_level * (i - 1) + 3, ...
        nline_Cmesh + (i - 1) * ncirc_line_per_level + 6, ...
        nline_per_level * (i - 1) + 5);      
    
%     fprintf(fo, 'CL_ext_%d = newcl;\n', nloop_per_level * (i - 1) + 4);
%     fprintf(fo, 'Curve Loop(CL_ext_%d) = {- %d, L_ext_%d, - %d, - L_ext_%d};\n', ...
%         nloop_per_level * (i - 1) + 4, ...
%         nline_airfoil + (i - 1) * n_ext_line_per_level + 5, ...
%         nline_per_level * (i - 1) + 4, ...
%         nline_Cmesh + (i - 1) * ncirc_line_per_level + 6, ...
%         nline_per_level * (i - 1) + 5);     
    
    fprintf(fo, 'CL_ext_%d = newcl;\n', nloop_per_level * (i - 1) + 5);
    fprintf(fo, 'Curve Loop(CL_ext_%d) = {%d, L_ext_%d, - %d, - L_ext_%d};\n', ...
        nloop_per_level * (i - 1) + 5, ...
        nline_airfoil + (i - 1) * n_ext_line_per_level + 6, ...
        nline_per_level * (i - 1) + 5, ...
        nline_Cmesh + (i - 1) * ncirc_line_per_level + 4, ...
        nline_per_level * (i - 1) + 7); 
    
%     fprintf(fo, 'CL_ext_%d = newcl;\n', nloop_per_level * (i - 1) + 6);
%     fprintf(fo, 'Curve Loop(CL_ext_%d) = {L_ext_%d, - %d, - L_ext_%d};\n', ...
%         nloop_per_level * (i - 1) + 6, ...
%         nline_per_level * (i - 1) + 6, ...
%         nline_Cmesh + (i - 1) * ncirc_line_per_level + 4, ...
%         nline_per_level * (i - 1) + 7);  
    
    fprintf(fo, 'CL_ext_%d = newcl;\n', nloop_per_level * (i - 1) + 7);
    fprintf(fo, 'Curve Loop(CL_ext_%d) = {%d, L_ext_%d, - %d, - L_ext_%d};\n', ...
        nloop_per_level * (i - 1) + 7, ...
        nline_airfoil + (i - 1) * n_ext_line_per_level + 4, ...
        nline_per_level * (i - 1) + 7, ...
        nline_Cmesh + (i - 1) * ncirc_line_per_level + 3, ...
        nline_per_level * (i - 1) + 8);     
    
    fprintf(fo, 'CL_ext_%d = newcl;\n', nloop_per_level * (i - 1) + 8);
    fprintf(fo, 'Curve Loop(CL_ext_%d) = {%d, L_ext_%d, - %d, - L_ext_%d};\n', ...
        nloop_per_level * (i - 1) + 8, ...
        nline_airfoil + (i - 1) * n_ext_line_per_level + 2, ...
        nline_per_level * (i - 1) + 8, ...
        nline_Cmesh + (i - 1) * ncirc_line_per_level + 2, ...
        nline_per_level * (i - 1) + 9);      
    
    fprintf(fo, 'CL_ext_%d = newcl;\n', nloop_per_level * (i - 1) + 9);
    fprintf(fo, 'Curve Loop(CL_ext_%d) = {- %d, L_ext_%d, - %d, - L_ext_%d};\n', ...
        nloop_per_level * (i - 1) + 9, ...
        nline_airfoil + (i - 1) * n_ext_line_per_level + 7, ...
        nline_per_level * (i - 1) + 9, ...
        nline_Cmesh + (i - 1) * ncirc_line_per_level + 1, ...
        nline_per_level * (i - 1) + 10);  
    
    fprintf(fo, 'CL_ext_%d = newcl;\n', nloop_per_level * (i - 1) + 10);
    fprintf(fo, 'Curve Loop(CL_ext_%d) = {%d, L_ext_%d, - %d, - L_ext_%d};\n', ...
        nloop_per_level * (i - 1) + 10, ...
        nline_airfoil + (i - 1) * n_ext_line_per_level + 8, ...
        nline_per_level * (i - 1) + 10, ...
        nline_Cmesh + (i - 1) * ncirc_line_per_level + 10, ...
        nline_per_level * (i - 1) + 1);     
    
    
    % Vertical faces
    
    fprintf(fo, 'CL_ext_%d = newcl;\n', nloop_per_level * (i - 1) + 11);
    fprintf(fo, 'Curve Loop(CL_ext_%d) = {- L_airfoil_%d, L_ext_%d, L_ext_%d, - L_ext_%d};\n', ...
        nloop_per_level * (i - 1) + 11, ...
        (i - 1) * nline_airfoil_per_level + 13, ...
        nline_per_level * (i - 1) + 1, ...
        nline_per_level * (i - 1) + 11, ...
        nline_per_level * i + 1);     
    
    fprintf(fo, 'CL_ext_%d = newcl;\n', nloop_per_level * (i - 1) + 12);
    fprintf(fo, 'Curve Loop(CL_ext_%d) = {- L_airfoil_%d, L_ext_%d, L_ext_%d, - L_ext_%d};\n', ...
        nloop_per_level * (i - 1) + 12, ...
        (i - 1) * nline_airfoil_per_level + 15, ...
        nline_per_level * (i - 1) + 2, ...
        nline_per_level * (i - 1) + 12, ...
        nline_per_level * i + 2);    
    
    fprintf(fo, 'CL_ext_%d = newcl;\n', nloop_per_level * (i - 1) + 13);
    fprintf(fo, 'Curve Loop(CL_ext_%d) = {- L_airfoil_%d, L_ext_%d, L_ext_%d, - L_ext_%d};\n', ...
        nloop_per_level * (i - 1) + 13, ...
        (i - 1) * nline_airfoil_per_level + 11, ...
        nline_per_level * (i - 1) + 3, ...
        nline_per_level * (i - 1) + 13, ...
        nline_per_level * i + 3);     
    
%     fprintf(fo, 'CL_ext_%d = newcl;\n', nloop_per_level * (i - 1) + 14);
%     fprintf(fo, 'Curve Loop(CL_ext_%d) = {- L_airfoil_%d, L_ext_%d, L_ext_%d, - L_ext_%d};\n', ...
%         nloop_per_level * (i - 1) + 14, ...
%         (i - 1) * nline_airfoil_per_level + 11, ...
%         nline_per_level * (i - 1) + 4, ...
%         nline_per_level * (i - 1) + 14, ...
%         nline_per_level * i + 4);      
    
    fprintf(fo, 'CL_ext_%d = newcl;\n', nloop_per_level * (i - 1) + 15);
    fprintf(fo, 'Curve Loop(CL_ext_%d) = {- L_airfoil_%d, L_ext_%d, L_ext_%d, - L_ext_%d};\n', ...
        nloop_per_level * (i - 1) + 15, ...
        (i - 1) * nline_airfoil_per_level + 17, ...
        nline_per_level * (i - 1) + 5, ...
        nline_per_level * (i - 1) + 15, ...
        nline_per_level * i + 5);     
    
%     fprintf(fo, 'CL_ext_%d = newcl;\n', nloop_per_level * (i - 1) + 16);
%     fprintf(fo, 'Curve Loop(CL_ext_%d) = {- L_airfoil_%d, L_ext_%d, L_ext_%d, - L_ext_%d};\n', ...
%         nloop_per_level * (i - 1) + 16, ...
%         (i - 1) * nline_airfoil_per_level + 12, ...
%         nline_per_level * (i - 1) + 6, ...
%         nline_per_level * (i - 1) + 16, ...
%         nline_per_level * i + 6);  
    
    fprintf(fo, 'CL_ext_%d = newcl;\n', nloop_per_level * (i - 1) + 17);
    fprintf(fo, 'Curve Loop(CL_ext_%d) = {- L_airfoil_%d, L_ext_%d, L_ext_%d, - L_ext_%d};\n', ...
        nloop_per_level * (i - 1) + 17, ...
        (i - 1) * nline_airfoil_per_level + 12, ...
        nline_per_level * (i - 1) + 7, ...
        nline_per_level * (i - 1) + 17, ...
        nline_per_level * i + 7);      
    
    fprintf(fo, 'CL_ext_%d = newcl;\n', nloop_per_level * (i - 1) + 18);
    fprintf(fo, 'Curve Loop(CL_ext_%d) = {- L_airfoil_%d, L_ext_%d, L_ext_%d, - L_ext_%d};\n', ...
        nloop_per_level * (i - 1) + 18, ...
        (i - 1) * nline_airfoil_per_level + 16, ...
        nline_per_level * (i - 1) + 8, ...
        nline_per_level * (i - 1) + 18, ...
        nline_per_level * i + 8);      
    
    fprintf(fo, 'CL_ext_%d = newcl;\n', nloop_per_level * (i - 1) + 19);
    fprintf(fo, 'Curve Loop(CL_ext_%d) = {- L_airfoil_%d, L_ext_%d, L_ext_%d, - L_ext_%d};\n', ...
        nloop_per_level * (i - 1) + 19, ...
        (i - 1) * nline_airfoil_per_level + 14, ...
        nline_per_level * (i - 1) + 9, ...
        nline_per_level * (i - 1) + 19, ...
        nline_per_level * i + 9);          
 
    fprintf(fo, 'CL_ext_%d = newcl;\n', nloop_per_level * (i - 1) + 20);
    fprintf(fo, 'Curve Loop(CL_ext_%d) = {- L_airfoil_%d, L_ext_%d, L_ext_%d, - L_ext_%d};\n', ...
        nloop_per_level * (i - 1) + 20, ...
        (i - 1) * nline_airfoil_per_level + 18, ...
        nline_per_level * (i - 1) + 10, ...
        nline_per_level * (i - 1) + 20, ...
        nline_per_level * i + 10);       
 
    
  % Exterior faces
  
    fprintf(fo, 'CL_ext_%d = newcl;\n', nloop_per_level * (i - 1) + 21);
    fprintf(fo, 'Curve Loop(CL_ext_%d) = {- L_ext_%d, - %d, L_ext_%d, %d};\n', ...
        nloop_per_level * (i - 1) + 21, ...
        nline_per_level * (i - 1) + 11, ...
        nline_Cmesh + (i - 1) * ncirc_line_per_level + 9, ...
        nline_per_level * (i - 1) + 12, ...
        nline_Cmesh + i * ncirc_line_per_level + 9);   
    
    fprintf(fo, 'CL_ext_%d = newcl;\n', nloop_per_level * (i - 1) + 22);
    fprintf(fo, 'Curve Loop(CL_ext_%d) = {- L_ext_%d, - %d, L_ext_%d, %d};\n', ...
        nloop_per_level * (i - 1) + 22, ...
        nline_per_level * (i - 1) + 12, ...
        nline_Cmesh + (i - 1) * ncirc_line_per_level + 8, ...
        nline_per_level * (i - 1) + 13, ...
        nline_Cmesh + i * ncirc_line_per_level + 8);      
    
%     fprintf(fo, 'CL_ext_%d = newcl;\n', nloop_per_level * (i - 1) + 23);
%     fprintf(fo, 'Curve Loop(CL_ext_%d) = {- L_ext_%d, - %d, L_ext_%d, %d};\n', ...
%         nloop_per_level * (i - 1) + 23, ...
%         nline_per_level * (i - 1) + 13, ...
%         nline_Cmesh + (i - 1) * ncirc_line_per_level + 7, ...
%         nline_per_level * (i - 1) + 14, ...
%         nline_Cmesh + i * ncirc_line_per_level + 7);     

    fprintf(fo, 'CL_ext_%d = newcl;\n', nloop_per_level * (i - 1) + 23);
    fprintf(fo, 'Curve Loop(CL_ext_%d) = {- L_ext_%d, - %d, L_ext_%d, %d};\n', ...
        nloop_per_level * (i - 1) + 23, ...
        nline_per_level * (i - 1) + 13, ...
        nline_Cmesh + (i - 1) * ncirc_line_per_level + 6, ...
        nline_per_level * (i - 1) + 15, ...
        nline_Cmesh + i * ncirc_line_per_level + 6);     

    
%     fprintf(fo, 'CL_ext_%d = newcl;\n', nloop_per_level * (i - 1) + 24);
%     fprintf(fo, 'Curve Loop(CL_ext_%d) = {- L_ext_%d, - %d, L_ext_%d, %d};\n', ...
%         nloop_per_level * (i - 1) + 24, ...
%         nline_per_level * (i - 1) + 14, ...
%         nline_Cmesh + (i - 1) * ncirc_line_per_level + 6, ...
%         nline_per_level * (i - 1) + 15, ...
%         nline_Cmesh + i * ncirc_line_per_level + 6);     
    
%     fprintf(fo, 'CL_ext_%d = newcl;\n', nloop_per_level * (i - 1) + 25);
%     fprintf(fo, 'Curve Loop(CL_ext_%d) = {- L_ext_%d, - %d, L_ext_%d, %d};\n', ...
%         nloop_per_level * (i - 1) + 25, ...
%         nline_per_level * (i - 1) + 15, ...
%         nline_Cmesh + (i - 1) * ncirc_line_per_level + 5, ...
%         nline_per_level * (i - 1) + 16, ...
%         nline_Cmesh + i * ncirc_line_per_level + 5);     
%     
%     fprintf(fo, 'CL_ext_%d = newcl;\n', nloop_per_level * (i - 1) + 26);
%     fprintf(fo, 'Curve Loop(CL_ext_%d) = {- L_ext_%d, - %d, L_ext_%d, %d};\n', ...
%         nloop_per_level * (i - 1) + 26, ...
%         nline_per_level * (i - 1) + 16, ...
%         nline_Cmesh + (i - 1) * ncirc_line_per_level + 4, ...
%         nline_per_level * (i - 1) + 17, ...
%         nline_Cmesh + i * ncirc_line_per_level + 4);    

    fprintf(fo, 'CL_ext_%d = newcl;\n', nloop_per_level * (i - 1) + 25);
    fprintf(fo, 'Curve Loop(CL_ext_%d) = {- L_ext_%d, - %d, L_ext_%d, %d};\n', ...
        nloop_per_level * (i - 1) + 25, ...
        nline_per_level * (i - 1) + 15, ...
        nline_Cmesh + (i - 1) * ncirc_line_per_level + 4, ...
        nline_per_level * (i - 1) + 17, ...
        nline_Cmesh + i * ncirc_line_per_level + 4); 
    
    fprintf(fo, 'CL_ext_%d = newcl;\n', nloop_per_level * (i - 1) + 27);
    fprintf(fo, 'Curve Loop(CL_ext_%d) = {- L_ext_%d, - %d, L_ext_%d, %d};\n', ...
        nloop_per_level * (i - 1) + 27, ...
        nline_per_level * (i - 1) + 17, ...
        nline_Cmesh + (i - 1) * ncirc_line_per_level + 3, ...
        nline_per_level * (i - 1) + 18, ...
        nline_Cmesh + i * ncirc_line_per_level + 3);     
    
    fprintf(fo, 'CL_ext_%d = newcl;\n', nloop_per_level * (i - 1) + 28);
    fprintf(fo, 'Curve Loop(CL_ext_%d) = {- L_ext_%d, - %d, L_ext_%d, %d};\n', ...
        nloop_per_level * (i - 1) + 28, ...
        nline_per_level * (i - 1) + 18, ...
        nline_Cmesh + (i - 1) * ncirc_line_per_level + 2, ...
        nline_per_level * (i - 1) + 19, ...
        nline_Cmesh + i * ncirc_line_per_level + 2);   
    
    fprintf(fo, 'CL_ext_%d = newcl;\n', nloop_per_level * (i - 1) + 29);
    fprintf(fo, 'Curve Loop(CL_ext_%d) = {- L_ext_%d, - %d, L_ext_%d, %d};\n', ...
        nloop_per_level * (i - 1) + 29, ...
        nline_per_level * (i - 1) + 19, ...
        nline_Cmesh + (i - 1) * ncirc_line_per_level + 1, ...
        nline_per_level * (i - 1) + 20, ...
        nline_Cmesh + i * ncirc_line_per_level + 1);     
    
    fprintf(fo, 'CL_ext_%d = newcl;\n', nloop_per_level * (i - 1) + 30);
    fprintf(fo, 'Curve Loop(CL_ext_%d) = {- L_ext_%d, - %d, L_ext_%d, %d};\n', ...
        nloop_per_level * (i - 1) + 30, ...
        nline_per_level * (i - 1) + 20, ...
        nline_Cmesh + (i - 1) * ncirc_line_per_level + 10, ...
        nline_per_level * (i - 1) + 11, ...
        nline_Cmesh + i * ncirc_line_per_level + 10);     
end

i = nradii;
fprintf(fo, 'CL_ext_%d = newcl;\n', nloop_per_level * (i - 1) + 1);
fprintf(fo, 'Curve Loop(CL_ext_%d) = {- %d, L_ext_%d, - %d, - L_ext_%d};\n', ...
    nloop_per_level * (i - 1) + 1, ...
    nline_airfoil + (i - 1) * n_ext_line_per_level + 1, ...
    nline_per_level * (i - 1) + 1, ...
    nline_Cmesh + (i - 1) * ncirc_line_per_level + 9, ...
    nline_per_level * (i - 1) + 2);      

fprintf(fo, 'CL_ext_%d = newcl;\n', nloop_per_level * (i - 1) + 2);
fprintf(fo, 'Curve Loop(CL_ext_%d) = {- %d, L_ext_%d, - %d, - L_ext_%d};\n', ...
    nloop_per_level * (i - 1) + 2, ...
    nline_airfoil + (i - 1) * n_ext_line_per_level + 3, ...
    nline_per_level * (i - 1) + 2, ...
    nline_Cmesh + (i - 1) * ncirc_line_per_level + 8, ...
    nline_per_level * (i - 1) + 3);    

fprintf(fo, 'CL_ext_%d = newcl;\n', nloop_per_level * (i - 1) + 3);
fprintf(fo, 'Curve Loop(CL_ext_%d) = {- %d, L_ext_%d, - %d, - L_ext_%d};\n', ...
    nloop_per_level * (i - 1) + 3, ...
    nline_airfoil + (i - 1) * n_ext_line_per_level + 5, ...
    nline_per_level * (i - 1) + 3, ...
    nline_Cmesh + (i - 1) * ncirc_line_per_level + 6, ...
    nline_per_level * (i - 1) + 5);      

%     fprintf(fo, 'CL_ext_%d = newcl;\n', nloop_per_level * (i - 1) + 4);
%     fprintf(fo, 'Curve Loop(CL_ext_%d) = {- %d, L_ext_%d, - %d, - L_ext_%d};\n', ...
%         nloop_per_level * (i - 1) + 4, ...
%         nline_airfoil + (i - 1) * n_ext_line_per_level + 5, ...
%         nline_per_level * (i - 1) + 4, ...
%         nline_Cmesh + (i - 1) * ncirc_line_per_level + 6, ...
%         nline_per_level * (i - 1) + 5);     

fprintf(fo, 'CL_ext_%d = newcl;\n', nloop_per_level * (i - 1) + 5);
fprintf(fo, 'Curve Loop(CL_ext_%d) = {%d, L_ext_%d, - %d, - L_ext_%d};\n', ...
    nloop_per_level * (i - 1) + 5, ...
    nline_airfoil + (i - 1) * n_ext_line_per_level + 6, ...
    nline_per_level * (i - 1) + 5, ...
    nline_Cmesh + (i - 1) * ncirc_line_per_level + 4, ...
    nline_per_level * (i - 1) + 7); 

%     fprintf(fo, 'CL_ext_%d = newcl;\n', nloop_per_level * (i - 1) + 6);
%     fprintf(fo, 'Curve Loop(CL_ext_%d) = {L_ext_%d, - %d, - L_ext_%d};\n', ...
%         nloop_per_level * (i - 1) + 6, ...
%         nline_per_level * (i - 1) + 6, ...
%         nline_Cmesh + (i - 1) * ncirc_line_per_level + 4, ...
%         nline_per_level * (i - 1) + 7);  

fprintf(fo, 'CL_ext_%d = newcl;\n', nloop_per_level * (i - 1) + 7);
fprintf(fo, 'Curve Loop(CL_ext_%d) = {%d, L_ext_%d, - %d, - L_ext_%d};\n', ...
    nloop_per_level * (i - 1) + 7, ...
    nline_airfoil + (i - 1) * n_ext_line_per_level + 4, ...
    nline_per_level * (i - 1) + 7, ...
    nline_Cmesh + (i - 1) * ncirc_line_per_level + 3, ...
    nline_per_level * (i - 1) + 8);     

fprintf(fo, 'CL_ext_%d = newcl;\n', nloop_per_level * (i - 1) + 8);
fprintf(fo, 'Curve Loop(CL_ext_%d) = {%d, L_ext_%d, - %d, - L_ext_%d};\n', ...
    nloop_per_level * (i - 1) + 8, ...
    nline_airfoil + (i - 1) * n_ext_line_per_level + 2, ...
    nline_per_level * (i - 1) + 8, ...
    nline_Cmesh + (i - 1) * ncirc_line_per_level + 2, ...
    nline_per_level * (i - 1) + 9);      

fprintf(fo, 'CL_ext_%d = newcl;\n', nloop_per_level * (i - 1) + 9);
fprintf(fo, 'Curve Loop(CL_ext_%d) = {- %d, L_ext_%d, - %d, - L_ext_%d};\n', ...
    nloop_per_level * (i - 1) + 9, ...
    nline_airfoil + (i - 1) * n_ext_line_per_level + 7, ...
    nline_per_level * (i - 1) + 9, ...
    nline_Cmesh + (i - 1) * ncirc_line_per_level + 1, ...
    nline_per_level * (i - 1) + 10);  

fprintf(fo, 'CL_ext_%d = newcl;\n', nloop_per_level * (i - 1) + 10);
fprintf(fo, 'Curve Loop(CL_ext_%d) = {%d, L_ext_%d, - %d, - L_ext_%d};\n', ...
    nloop_per_level * (i - 1) + 10, ...
    nline_airfoil + (i - 1) * n_ext_line_per_level + 8, ...
    nline_per_level * (i - 1) + 10, ...
    nline_Cmesh + (i - 1) * ncirc_line_per_level + 10, ...
    nline_per_level * (i - 1) + 1);         



%% Declare the surfaces

fprintf(fo, '\n' );
for i = 1 : nradii - 1
    for j = 1 : nloop_per_level
%     for j = 1 : 20
        if(j ~= 4 && j~= 6 && j~= 14 && j~= 16 && j~= 24 && j~= 26)
            fprintf(fo, 'S_ext_%d = news;\n', (i - 1) * nloop_per_level + j); 
            fprintf(fo, 'Surface(S_ext_%d) = {CL_ext_%d};\n', ...
                (i - 1) * nloop_per_level + j, (i - 1) * nloop_per_level + j);
        end
    end
end
i = nradii;
for j = 1 : 10
    if(j ~= 4 && j~= 6 && j~= 14 && j~= 16 && j~= 24 && j~= 26)
        fprintf(fo, 'S_ext_%d = news;\n', (i - 1) * nloop_per_level + j); 
        fprintf(fo, 'Surface(S_ext_%d) = {CL_ext_%d};\n', ...
            (i - 1) * nloop_per_level + j, (i - 1) * nloop_per_level + j);
    end
end


%% Make the surfaces transfinite (careful, not the ones in the corner)

fprintf(fo, '\n' );
for i = 1 : nradii - 1
    for j = 1 : nloop_per_level
        if(j ~= 4 && j~= 6 && j~= 14 && j~= 16 && j~= 24 && j~= 26)
            fprintf(fo, 'Transfinite Surface {S_ext_%d};\n', (i - 1) * nloop_per_level + j);      
        end
    end
end
i = nradii;
for j = 1 : 10
    if(j ~= 4 && j~= 6 && j~= 14 && j~= 16 && j~= 24 && j~= 26)
        fprintf(fo, 'Transfinite Surface {S_ext_%d};\n', (i - 1) * nloop_per_level + j);  
    end
end


%% Recombine the surfaces (but not the ones in the corner)

fprintf(fo, '\n' );
fprintf(fo, 'Recombine Surface {');
for i = 1 : nradii - 1
    for j = 1 : nloop_per_level
        if(j ~= 4 && j~= 6 && j~= 14 && j~= 16 && j~= 24 && j~= 26)
            fprintf(fo, 'S_ext_%d, ', (i - 1) * nloop_per_level + j);
        end
    end
end
i = nradii;
for j = 1 : 9
    if(j ~= 4 && j~= 6 && j~= 14 && j~= 16 && j~= 24 && j~= 26)
        fprintf(fo, 'S_ext_%d, ', (i - 1) * nloop_per_level + j);
    end
end
fprintf(fo, 'S_ext_%d', (i - 1) * nloop_per_level + 10);
fprintf(fo, '};\n' );


%% Declare the surface loops

nsloop_per_level = 10;
nsurf_airfoil_per_level = 24;
fprintf(fo, '\n' );


for i = 1 : nradii - 1
    fprintf(fo, 'SL_ext_%d = newsl;\n', (i - 1) * nsloop_per_level + 1);
    fprintf(fo, 'Surface Loop(SL_ext_%d) = {S_ext_%d, S_ext_%d, S_ext_%d, S_ext_%d, S_ext_%d, S_airfoil_%d};\n', ...
        (i - 1) * nsloop_per_level + 1, ...
        nloop_per_level * (i - 1) + 1, ...
        nloop_per_level * (i - 1) + 11, ...
        nloop_per_level * (i - 1) + 12 , ...
        nloop_per_level * i + 1, ...
        nloop_per_level * (i - 1) + 21, ...
        nsurf_airfoil_per_level * (i - 1) + 13); 


    fprintf(fo, 'SL_ext_%d = newsl;\n', (i - 1) * nsloop_per_level + 2);
    fprintf(fo, 'Surface Loop(SL_ext_%d) = {S_ext_%d, S_ext_%d, S_ext_%d, S_ext_%d, S_ext_%d, S_airfoil_%d};\n', ...
        (i - 1) * nsloop_per_level + 2, ...
        nloop_per_level * (i - 1) + 2, ...
        nloop_per_level * (i - 1) + 12, ...
        nloop_per_level * (i - 1) + 13 , ...
        nloop_per_level * i + 2, ...
        nloop_per_level * (i - 1) + 22, ...
        nsurf_airfoil_per_level * (i - 1) + 14); 



%     fprintf(fo, 'SL_ext_%d = newsl;\n', (i - 1) * nsloop_per_level + 3);
%     fprintf(fo, 'Surface Loop(SL_ext_%d) = {S_ext_%d, S_ext_%d, S_ext_%d, S_ext_%d, S_ext_%d};\n', ...
%         (i - 1) * nsloop_per_level + 3, ...
%         nloop_per_level * (i - 1) + 3, ...
%         nloop_per_level * (i - 1) + 13, ...
%         nloop_per_level * (i - 1) + 14, ...
%         nloop_per_level * i + 3, ...
%         nloop_per_level * (i - 1) + 23);
    
%     fprintf(fo, 'SL_ext_%d = newsl;\n', (i - 1) * nsloop_per_level + 4);
%     fprintf(fo, 'Surface Loop(SL_ext_%d) = {S_ext_%d, S_ext_%d, S_ext_%d, S_ext_%d, S_ext_%d, S_airfoil_%d};\n', ...
%         (i - 1) * nsloop_per_level + 4, ...
%         nloop_per_level * (i - 1) + 4, ...
%         nloop_per_level * (i - 1) + 14, ...
%         nloop_per_level * (i - 1) + 15, ...
%         nloop_per_level * i + 4, ...
%         nloop_per_level * (i - 1) + 24, ...
%         nsurf_airfoil_per_level * (i - 1) + 15);
    

    fprintf(fo, 'SL_ext_%d = newsl;\n', (i - 1) * nsloop_per_level + 4);
    fprintf(fo, 'Surface Loop(SL_ext_%d) = {S_ext_%d, S_ext_%d, S_ext_%d, S_ext_%d, S_ext_%d, S_airfoil_%d};\n', ...
        (i - 1) * nsloop_per_level + 4, ...
        nloop_per_level * (i - 1) + 3, ...
        nloop_per_level * (i - 1) + 13, ...
        nloop_per_level * (i - 1) + 15, ...
        nloop_per_level * i + 3, ...
        nloop_per_level * (i - 1) + 23, ...
        nsurf_airfoil_per_level * (i - 1) + 15);

%     fprintf(fo, 'SL_ext_%d = newsl;\n', (i - 1) * nsloop_per_level + 5);
%     fprintf(fo, 'Surface Loop(SL_ext_%d) = {S_ext_%d, S_ext_%d, S_ext_%d, S_ext_%d, S_ext_%d, S_airfoil_%d};\n', ...
%         (i - 1) * nsloop_per_level + 5, ...
%         nloop_per_level * (i - 1) + 5, ...
%         nloop_per_level * (i - 1) + 15, ...
%         nloop_per_level * (i - 1) + 16, ...
%         nloop_per_level * i + 5, ...
%         nloop_per_level * (i - 1) + 25, ...
%         nsurf_airfoil_per_level * (i - 1) + 16);    

    fprintf(fo, 'SL_ext_%d = newsl;\n', (i - 1) * nsloop_per_level + 5);
    fprintf(fo, 'Surface Loop(SL_ext_%d) = {S_ext_%d, S_ext_%d, S_ext_%d, S_ext_%d, S_ext_%d, S_airfoil_%d};\n', ...
        (i - 1) * nsloop_per_level + 5, ...
        nloop_per_level * (i - 1) + 5, ...
        nloop_per_level * (i - 1) + 15, ...
        nloop_per_level * (i - 1) + 17, ...
        nloop_per_level * i + 5, ...
        nloop_per_level * (i - 1) + 25, ...
        nsurf_airfoil_per_level * (i - 1) + 16);    
    

%     fprintf(fo, 'SL_ext_%d = newsl;\n', (i - 1) * nsloop_per_level + 6);
%     fprintf(fo, 'Surface Loop(SL_ext_%d) = {S_ext_%d, S_ext_%d, S_ext_%d, S_ext_%d, S_ext_%d};\n', ...
%         (i - 1) * nsloop_per_level + 6, ...
%         nloop_per_level * (i - 1) + 6, ...
%         nloop_per_level * (i - 1) + 16, ...
%         nloop_per_level * (i - 1) + 17, ...
%         nloop_per_level * i + 6, ...
%         nloop_per_level * (i - 1) + 26);
    
    fprintf(fo, 'SL_ext_%d = newsl;\n', (i - 1) * nsloop_per_level + 7);
    fprintf(fo, 'Surface Loop(SL_ext_%d) = {S_ext_%d, S_ext_%d, S_ext_%d, S_ext_%d, S_ext_%d, S_airfoil_%d};\n', ...
        (i - 1) * nsloop_per_level + 7, ...
        nloop_per_level * (i - 1) + 7, ...
        nloop_per_level * (i - 1) + 17, ...
        nloop_per_level * (i - 1) + 18, ...
        nloop_per_level * i + 7, ...
        nloop_per_level * (i - 1) + 27, ...
        nsurf_airfoil_per_level * (i - 1) + 17);    
    
    fprintf(fo, 'SL_ext_%d = newsl;\n', (i - 1) * nsloop_per_level + 8);
    fprintf(fo, 'Surface Loop(SL_ext_%d) = {S_ext_%d, S_ext_%d, S_ext_%d, S_ext_%d, S_ext_%d, S_airfoil_%d};\n', ...
        (i - 1) * nsloop_per_level + 8, ...
        nloop_per_level * (i - 1) + 8, ...
        nloop_per_level * (i - 1) + 18, ...
        nloop_per_level * (i - 1) + 19, ...
        nloop_per_level * i + 8, ...
        nloop_per_level * (i - 1) + 28, ...
        nsurf_airfoil_per_level * (i - 1) + 18);  
    
    fprintf(fo, 'SL_ext_%d = newsl;\n', (i - 1) * nsloop_per_level + 9);
    fprintf(fo, 'Surface Loop(SL_ext_%d) = {S_ext_%d, S_ext_%d, S_ext_%d, S_ext_%d, S_ext_%d, S_airfoil_%d};\n', ...
        (i - 1) * nsloop_per_level + 9, ...
        nloop_per_level * (i - 1) + 9, ...
        nloop_per_level * (i - 1) + 19, ...
        nloop_per_level * (i - 1) + 20, ...
        nloop_per_level * i + 9, ...
        nloop_per_level * (i - 1) + 29, ...
        nsurf_airfoil_per_level * (i - 1) + 19);  
    
    fprintf(fo, 'SL_ext_%d = newsl;\n', (i - 1) * nsloop_per_level + 10);
    fprintf(fo, 'Surface Loop(SL_ext_%d) = {S_ext_%d, S_ext_%d, S_ext_%d, S_ext_%d, S_ext_%d, S_airfoil_%d};\n', ...
        (i - 1) * nsloop_per_level + 10, ...
        nloop_per_level * (i - 1) + 10, ...
        nloop_per_level * (i - 1) + 20, ...
        nloop_per_level * (i - 1) + 11, ...
        nloop_per_level * i + 10, ...
        nloop_per_level * (i - 1) + 30, ...
        nsurf_airfoil_per_level * (i - 1) + 20);      
end

%% Declare the volumes

nvolume_per_level = 10;

fprintf(fo, '\n' );
for i = 1 : nradii - 1
    for j = 1 : nvolume_per_level
        if(j ~= 3 && j~= 6)
            fprintf(fo, 'V_ext_%d = newv;\n', (i - 1) * nvolume_per_level + j);
            fprintf(fo, 'Volume(V_ext_%d) = {SL_ext_%d};\n', (i - 1) * nvolume_per_level + j, (i - 1) * nsloop_per_level + j);
        end
    end
end

%% Make the volumes transfinite (careful, not the ones in the corner)

fprintf(fo, '\n' );
for i = 1 : nradii - 1
    for j = 1 : nvolume_per_level
        if(j ~= 3 && j~= 6)
            fprintf(fo, 'Transfinite Volume {V_ext_%d};\n', (i - 1) * nvolume_per_level + j);
        end
    end
end


end

