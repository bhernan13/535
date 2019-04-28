function output_properties = material_properties(material)

% This function has all of the researched properties for potential pebble
% bed materials. It will output the properties for the selected material.

%% all material properties
properties = cell(6,4);
properties{1,1} = 'Material'; properties{1,2} = 'T [K];k [W/m-k]';
properties{1,3} = 'T [K]; Cp [J/kg-K]'; properties{1,4} = 'Density [lbm/ft3]';
properties{2,1} = 'SS304'; properties{3,1} = 'SS316';
properties{4,1} = 'CU201'; properties{5,1} = 'AL2017';
properties{6,1} = 'BRASS260';
% Thermal Conductivity 
properties{2,2} = [ 293, 373, 473, 573, 673, 773;
                     16.2, 16.2, 17.5, 18.8, 20.1, 21.4];
properties{3,2} = [293, 373, 473, 573, 673, 773;
                   14, 14.9, 16, 17.3, 18.6, 19.9];
properties{4,2} = [250, 300, 350, 400, 500, 600, 800;
                   406, 401, 396, 393, 386, 379, 366];
properties{5,2} = [250, 300, 400, 500, 600, 700, 800
                   134, 135, 136, 135, 132, 132, 125];
properties{6,2} = [293, 350, 400, 500, 600, 700, 800
                   120, 119, 118, 116, 114, 112, 11];
% Specific Heat - Put everything into cell arrays
properties{2,3} = @(T) 6.683 + 0.04906*T + 80.74*log(T); %{'A + BT + Cln(T)','Equation';'A', 6.683;'B', 0.04906;'C', 80.74};

properties{3,3} = [293, 363, 473, 593, 703, 813, 923, 1033, 1143;
                   452, 486, 528, 548, 565, 573, 586, 615, 649];
properties{4,3} = [300, 400, 600, 800;
                   385, 397, 417, 433];
properties{5,3} = [298, 373, 473, 573, 673, 773, 811;
                   850, 900, 950, 970, 1000, 1080, 1100];
properties{6,3} = [300, 400, 600;
                   380, 395, 425];
% Density
properties{2,4} = 499.392;
properties{3,4} = 499.392;
properties{4,4} = 558.144;
properties{5,4} = 174.528;
properties{6,4} = 532.224;

%% pick desired properties
materialSelect = strcmpi(properties(:,1),material);
if materialSelect(:) == 0
    error('Material data is not entered for chosen material. Choose a material with existing data or enter new material data.');
end
indexMaterial = find(materialSelect ~= 0);
Data_k = properties{indexMaterial,2};
Data_Cp = properties{indexMaterial,3};
Data_rho = properties{indexMaterial,4};

output_properties = cell(4,1);
output_properties{1} = material;
output_properties{2} = Data_k;
output_properties{3} = Data_Cp;
output_properties{4} = Data_rho;
end