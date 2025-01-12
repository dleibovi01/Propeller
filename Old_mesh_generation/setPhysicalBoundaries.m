function setPhysicalBoundaries(fo, radius, level_min, level_max)


nradii = length(radius);

%% Inlet 
fprintf(fo, '\n' );
fprintf(fo, '\n' );


nloop_per_level = 30;
fprintf(fo, 'Physical Surface("inlet") ={');
 
 for i = 1 : nradii - 2
     fprintf(fo, 'S_ext_%d, ',  nloop_per_level * (i - 1) + 25);
     fprintf(fo, 'S_ext_%d, ',  nloop_per_level * (i - 1) + 27);
     fprintf(fo, 'S_ext_%d, ',  nloop_per_level * (i - 1) + 28);
     fprintf(fo, 'S_ext_%d, ',  nloop_per_level * (i - 1) + 29);
 end
 i = nradii - 1; 
 fprintf(fo, 'S_ext_%d, ',  nloop_per_level * (i - 1) + 25);
 fprintf(fo, 'S_ext_%d, ',  nloop_per_level * (i - 1) + 27);
 fprintf(fo, 'S_ext_%d, ',  nloop_per_level * (i - 1) + 28);
 fprintf(fo, 'S_ext_%d};\n',  nloop_per_level * (i - 1) + 29);


%% Outlet

fprintf(fo, 'Physical Surface("outlet") ={');
 
 for i = 1 : nradii - 2
     fprintf(fo, 'S_ext_%d, ',  nloop_per_level * (i - 1) + 21);
     fprintf(fo, 'S_ext_%d, ',  nloop_per_level * (i - 1) + 22);
     fprintf(fo, 'S_ext_%d, ',  nloop_per_level * (i - 1) + 23);
     fprintf(fo, 'S_ext_%d, ',  nloop_per_level * (i - 1) + 30);
 end
 i = nradii - 1; 
 fprintf(fo, 'S_ext_%d, ',  nloop_per_level * (i - 1) + 21);
 fprintf(fo, 'S_ext_%d, ',  nloop_per_level * (i - 1) + 22);
 fprintf(fo, 'S_ext_%d, ',  nloop_per_level * (i - 1) + 23);
 fprintf(fo, 'S_ext_%d};\n',  nloop_per_level * (i - 1) + 30);



%% Bottom

fprintf(fo, '\n' );
fprintf(fo, 'Physical Surface("bottom") ={');
% 12 surfaces
fprintf(fo, 'S_ext_%d, S_ext_%d, S_ext_%d, S_ext_%d, S_ext_%d, S_ext_%d, S_ext_%d, S_ext_%d, S_airfoil_%d, S_airfoil_%d, S_airfoil_%d, S_airfoil_%d, S_airfoil_%d, S_airfoil_%d, S_blade_%d, S_blade_%d};\n', ...
    1, 2, 3, 5, 7, 8, 9, 10, 1, 2, 3, 4, 5, 6, 1, 2);


%% Top

nloop_per_level_Cmesh = 24;
fprintf(fo, '\n' );
fprintf(fo, 'Physical Surface("top") ={');
% 12 surfaces
fprintf(fo, 'S_ext_%d, S_ext_%d, S_ext_%d, S_ext_%d, S_ext_%d, S_ext_%d, S_ext_%d, S_ext_%d, S_airfoil_%d, S_airfoil_%d, S_airfoil_%d, S_airfoil_%d, S_airfoil_%d, S_airfoil_%d, S_blade_%d, S_blade_%d};\n', ...
    nloop_per_level * (nradii - 1) + 1, nloop_per_level * (nradii - 1) + 2, ...
    nloop_per_level * (nradii - 1) + 3, nloop_per_level * (nradii - 1) + 5, ...
    nloop_per_level * (nradii - 1) + 7, nloop_per_level * (nradii - 1) + 8, ...
    nloop_per_level * (nradii - 1) + 9, nloop_per_level * (nradii - 1) + 10, ...
    nloop_per_level_Cmesh * (nradii - 1) + 1, nloop_per_level_Cmesh * (nradii - 1) + 2, ...
    nloop_per_level_Cmesh * (nradii - 1) + 3, nloop_per_level_Cmesh * (nradii - 1) + 4, ...
    nloop_per_level_Cmesh * (nradii - 1) + 5, nloop_per_level_Cmesh * (nradii - 1) + 6, ...
    3 * (nradii - 1) + 1, 3 * (nradii - 1) + 2);
 
 

%% Blade

% fprintf(fo, '\n' );
% fprintf(fo, 'Physical Surface("blade") ={');
% nloop_per_level = 24;
% %  4 faces between each radii
% for i = level_min : level_max - 1
%     fprintf(fo, 'S_airfoil_%d, S_airfoil_%d, S_airfoil_%d, S_airfoil_%d,', ...
%         nloop_per_level * (i - 1) + 21, nloop_per_level * (i - 1) + 22, ...
%         nloop_per_level * (i - 1) + 23, nloop_per_level * (i - 1) + 24);
% end
% % then 4 faces, 2 a the bottom and 2 at the top
% fprintf(fo, 'S_blade_%d, S_blade_%d, S_blade_%d, S_blade_%d};\n', ...
%     3 * (level_min - 1) + 1, 3 * (level_min - 1) + 2, ...
%     3 * (level_max - 1) + 1, 3 * (level_max - 1) + 2);    


fprintf(fo, '\n' );
fprintf(fo, 'Physical Surface("blade_inlet") = {');
nloop_per_level = 24;
%  4 faces between each radii
for i = level_min : level_max - 2
    fprintf(fo, 'S_airfoil_%d, S_airfoil_%d,', ...
        nloop_per_level * (i - 1) + 21, nloop_per_level * (i - 1) + 24);
end
fprintf(fo, 'S_airfoil_%d, S_airfoil_%d};\n', ...
    nloop_per_level * (level_max - 2) + 21, nloop_per_level * (level_max - 2) + 24);


fprintf(fo, '\n' );
fprintf(fo, 'Physical Surface("blade_outlet") = {');
nloop_per_level = 24;
%  4 faces between each radii
for i = level_min : level_max - 2
    fprintf(fo, 'S_airfoil_%d, S_airfoil_%d,', ...
        nloop_per_level * (i - 1) + 22, nloop_per_level * (i - 1) + 23);
end
fprintf(fo, 'S_airfoil_%d, S_airfoil_%d};\n', ...
    nloop_per_level * (level_max - 2) + 22, nloop_per_level * (level_max - 2) + 23);


% then 2 faces at the bottom
fprintf(fo, '\n' );
fprintf(fo, 'Physical Surface("blade_bottom") = {');
fprintf(fo, 'S_blade_%d, S_blade_%d};\n', 3 * (level_min - 1) + 1, 3 * (level_min - 1) + 2);   
 


%  and 2 faces at the top
fprintf(fo, '\n' );
fprintf(fo, 'Physical Surface("blade_top") = {');
fprintf(fo, 'S_blade_%d, S_blade_%d};\n',  3 * (level_max - 1) + 1, 3 * (level_max - 1) + 2);
   


%% Physical volume

nvolume_per_level = 10;
fprintf(fo, 'Physical Volume("domain") ={');
% all the exterior volumes
for i = 1 : nradii - 1
    for j = 1 : nvolume_per_level
        if(j ~= 3 && j~= 6)
            fprintf(fo, 'V_ext_%d, ', (i - 1) * nvolume_per_level + j);
        end
    end
end

% all the Cmesh volumes
nvolume_per_level = 6;
for i = 1 : nradii - 1
    for j = 1 : nvolume_per_level
        fprintf(fo, 'V_airfoil_%d, ', (i - 1) * nvolume_per_level + j);
    end
end

% all the interior volumes below level_min and above level_max
for i = 1 : nradii - 2
    if(i < level_min || i >= level_max)
        for j = 1 : 2
            fprintf(fo, 'V_blade_%d, ', 2 * (i - 1) + j);
        end
    end
end
fprintf(fo, 'V_blade_%d, ', 2 * (nradii - 2) + 1);
fprintf(fo, 'V_blade_%d};\n ', 2 * (nradii - 2) + 2);



end