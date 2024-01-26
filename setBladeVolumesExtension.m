

function [npoint, nline, nloop, nsurface] = setBladeVolumesExtension(fo, Ni, ...
    npoint, nline, nloop, nsurface, nradii, nres, ...
    n_airfoil_lateral, radius, rmin, rmax)



% NACA = [0 0 1 2];           % NACA 4-digit designation as a row vector;

sizeNACA = 5e-3;            % resolution at wall
sizeWake = 0.25; % 5e-2;            % resolution at wake
sizeFar  = 1;               % resolution at far field
shrink = 0.7;


n_vertical = nres * (60*shrink+25);
r_vertical=1.015;
% n_wake = nres * 110;

r_wake_airfoil = 1.005;

% n_vert = nres * 70;

% n_inlet = nres * 40;
r_inlet = 0.3;


Nsplit = round(0.24 * Ni);

ext_points_per_level = 18;

tol = 1e-8;


% Create the line to cut the airfoils in two

for i = 1 : nradii
    
    if(radius(i) < rmin + tol || radius(i) > rmax - tol)
        % newcl
        nloop_per_level = 24;
        n_ext_line_per_level = 8;
        nline_airfoil_per_level = 6;
        fprintf(fo, '\n' );



        % Horizontal faces    

        fprintf(fo, 'L_int_%d = newl;\n', i);
        fprintf(fo, 'Line(L_int_%d) = {%d, %d};\n', ...
            i, ...
            (i - 1) * 2 * (Ni - 1) + Nsplit, ...
            (i - 1) * 2 * (Ni - 1) + 2 * Ni - Nsplit);    
        
        fprintf(fo, '\n' );
        fprintf(fo, 'Transfinite Curve {L_int_%d} = %d Using Progression %f;\n', ...
            i, n_airfoil_lateral, 1);        
    end
    
end



% Create the curve loops

for i = 1 : nradii
    
    % newcl
    nloop_per_level = 24;
    n_ext_line_per_level = 8;
    nline_per_level = 18;
    nline_airfoil_per_level = 6; 
    
    if(radius(i) < rmin + tol || radius(i) > rmax - tol)

        fprintf(fo, '\n' );



        % Horizontal faces    
     

        fprintf(fo, 'CL_blade_%d = newcl;\n', 3 * (i - 1) + 1);
        fprintf(fo, 'Curve Loop(CL_blade_%d) = {%d, %d, L_int_%d};\n', ...
            3 * (i - 1) + 1, ...
            nline_airfoil_per_level * (i - 1) + 3, ...
            nline_airfoil_per_level * (i - 1) + 6, ...
            i);    
        
        fprintf(fo, 'CL_blade_%d = newcl;\n', 3 * (i - 1) + 2);
        fprintf(fo, 'Curve Loop(CL_blade_%d) = {- %d, - %d, L_int_%d};\n', ...
            3 * (i - 1) + 2, ...
            nline_airfoil_per_level * (i - 1) + 5, ...
            nline_airfoil_per_level * (i - 1) + 4, ...
            i);         
    end
    
    % Vertical faces  
    if(i < nradii)
        if(radius(i + 1) < rmin + tol || radius(i) > rmax - tol)
            fprintf(fo, 'CL_blade_%d = newcl;\n', 3 * (i - 1) + 3);
            fprintf(fo, 'Curve Loop(CL_blade_%d) = {L_int_%d, L_airfoil_%d, -L_int_%d, - L_airfoil_%d};\n', ...
                3 * (i - 1) + 3, ...
                i, ...
                nline_per_level * (i - 1) + 9, ...
                i + 1, ...
                nline_per_level * (i - 1) + 7);        
        end
    end

    
end



% Create the surfaces
% news
% fprintf(fo, 'Plane Surface(%d) = {%d};\n', nsurface + 1, nloop + 1);
fprintf(fo, '\n' );
for i = 1 : nradii
    if(radius(i) < rmin + tol || radius(i) > rmax - tol)
        fprintf(fo, 'S_blade_%d = news;\n', 3 * (i - 1) + 1); 
        fprintf(fo, 'Plane Surface(S_blade_%d) = {CL_blade_%d};\n', 3 * (i - 1) + 1, 3 * (i - 1) + 1);
        
        fprintf(fo, 'S_blade_%d = news;\n', 3 * (i - 1) + 2); 
        fprintf(fo, 'Plane Surface(S_blade_%d) = {CL_blade_%d};\n', 3 * (i - 1) + 2, 3 * (i - 1) + 2);        
    end
    
    if(i < nradii)
        if(radius(i + 1) < rmin + tol || radius(i) > rmax - tol)
            fprintf(fo, 'S_blade_%d = news;\n', 3 * (i - 1) + 3); 
            fprintf(fo, 'Plane Surface(S_blade_%d) = {CL_blade_%d};\n', 3 * (i - 1) + 3, 3 * (i - 1) + 3);        
        end
    end
end

% Make the surfaces transfinite

fprintf(fo, '\n' );
for i = 1 : nradii
    if(radius(i) < rmin + tol || radius(i) > rmax - tol)
        fprintf(fo, 'Transfinite Surface {S_blade_%d};\n', 3 * (i - 1) + 1);        
        fprintf(fo, 'Transfinite Surface {S_blade_%d};\n', 3 * (i - 1) + 2); 
    end
    if(i < nradii)
        if(radius(i + 1) < rmin + tol || radius(i) > rmax - tol)
            fprintf(fo, 'Transfinite Surface {S_blade_%d};\n', 3 * (i - 1) + 3);       
        end
    end    
end

% Recombine the surfaces

fprintf(fo, '\n' );

for i = 1 : nradii
    if(radius(i) < rmin + tol || radius(i) > rmax - tol)
        fprintf(fo, 'Recombine Surface {S_blade_%d};\n', 3 * (i - 1) + 1);
        fprintf(fo, 'Recombine Surface {S_blade_%d};\n', 3 * (i - 1) + 2);
    end
    if(i < nradii)
        if(radius(i + 1) < rmin + tol || radius(i) > rmax - tol)
            fprintf(fo, 'Recombine Surface {S_blade_%d};\n', 3 * (i - 1) + 3);    
        end
    end     
end


% 
% Create the surface loops
% newsl
nsloop_per_level = 6;
fprintf(fo, '\n' );
for i = 1 : nradii - 1
    if(radius(i + 1) < rmin + tol || radius(i) > rmax - tol)
        fprintf(fo, 'SL_blade_%d = newsl;\n', 2 * (i - 1) + 1);
        fprintf(fo, 'Surface Loop(SL_blade_%d) = {S_blade_%d, S_blade_%d, S_blade_%d, S_airfoil_%d, S_airfoil_%d};\n', ...
            2 * (i - 1) + 1, ...
            3 * (i - 1) + 1, ...
            3 * (i - 1) + 3, ...
            3 * i + 1, ...
            nloop_per_level * (i - 1) + 23, ...
            nloop_per_level * (i - 1) + 24);   
        
        fprintf(fo, 'SL_blade_%d = newsl;\n', 2 * (i - 1) + 2);
        fprintf(fo, 'Surface Loop(SL_blade_%d) = {S_blade_%d, S_blade_%d, S_blade_%d, S_airfoil_%d, S_airfoil_%d};\n', ...
            2 * (i - 1) + 2, ...
            3 * (i - 1) + 2, ...
            3 * (i - 1) + 3, ...
            3 * i + 2, ...
            nloop_per_level * (i - 1) + 21, ...
            nloop_per_level * (i - 1) + 22);          
    end
end


% fprintf(fo, 'V_blade_%d = newv;\n', 1);
% fprintf(fo, 'Volume(V_blade_%d) = {SL_blade_%d};\n', 1, 1);
% fprintf(fo, 'V_blade_%d = newv;\n', 2);
% fprintf(fo, 'Volume(V_blade_%d) = {SL_blade_%d};\n', 2, 2);





% Create the volumes
% newv
nvolume_per_level = 6;
fprintf(fo, '\n' );
% for i = 1 : 1
for i = 1 : nradii - 1
    if(radius(i + 1) < rmin + tol || radius(i) > rmax - tol)
        fprintf(fo, 'V_blade_%d = newv;\n', 2 * (i - 1) + 1);
        fprintf(fo, 'Volume(V_blade_%d) = {SL_blade_%d};\n', 2 * (i - 1) + 1, 2 * (i - 1) + 1);
        
        fprintf(fo, 'V_blade_%d = newv;\n', 2 * (i - 1) + 2);
        fprintf(fo, 'Volume(V_blade_%d) = {SL_blade_%d};\n', 2 * (i - 1) + 2, 2 * (i - 1) + 2);        
    end
end



% Create transfinite volumes

fprintf(fo, '\n' );
for i = 1 : nradii - 1    
    if(radius(i + 1) < rmin + tol || radius(i) > rmax - tol)
        fprintf(fo, 'Transfinite Volume {V_blade_%d};\n', 2 * (i - 1) + 1);
        fprintf(fo, 'Transfinite Volume {V_blade_%d};\n', 2 * (i - 1) + 2);       
    end    
end


end

