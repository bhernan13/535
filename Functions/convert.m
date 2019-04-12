classdef convert
    %% convert is a conversion class object which contains static methods.
    % The static methods contained in this class are for converting units
    % between imperial (English) and SI unit systems or vice versa.
    %
    % The following methods are included:
    %
    % - Length
    % - Area
    % - Volume
    % - Pressure
    % - Temperature
    % - Conduction
    % - Convection
    % - Force
    % - Density
    % - Flow
    % - Heat Capacity
    %
    % All methods have follow the scheme below.
    %
    % Inputs:
    % * value   -> Value to be converted
    % * inUnit  -> Unit of Input
    % * outUnit -> Unit of Output
    %
    % This function then checks user input for validity before
    % finding the ratio of the given inUnit to the standard unit of
    % each system. This function defines feet (ft) and meters (m)
    % as the standard units of measure for length in the imperial
    % and SI unit systems respectively.
    %
    % If a conversion between imperial and SI units are specified,
    % the conversion factor is defined as the ratio between the two
    % standard units in the direction of conversion.
    %
    % Vector inputs are accepted, but units in must be consistent and units
    % out must be consistent. Input vector in ft2, request answer in mm2.
    % Cannot input [ft2,m2] and get [mm2, mm2]. Can also not input        
    % [ft2,ft2] and get [m2,mm2]. See Example 3.
    %
    % To learn more about any method listed above, type help
    % convert.method, replacing method with a method listed above, to find
    % out more.
    %
    % Examples:
    %
    % [outValue, outUnit] = convert.length(1,'ft','m')
    %   Returns the value of feet to meters
    %
    % [outValue, outUnit] = convert.convection(value,'BTU/ft2-hr-R','W/m2-C')
    %   Returns the value of BTU/ft2-hr-R in W/m2-C
    %
    % [outValue, outUnit] = convert.area([1.5,2.3],'ft2','mm2')
    %   Returns the area of [input1,input2,...] ft2 [output1,output2,...] mm2. 
    %
    
    methods(Static)
        function [outVal, outUnit] = length(value,inUnit,outUnit)
            %% [outVal, outUnit] = length(value,inUnit,outUnit)
            % This function converts lengths of measurements from the
            % imperial (English) system to the SI system. 
            %
            %
            % ENG -> SI: 1 ft = 0.3048 m = conv
            % SI -> ENG: 1 m = 3.2808 ft = conv
            % SI -> SI or ENG -> ENG:  1 = conv
            % 
            % Supported unit inputs are:
            %
            % ENG: 'feet','ft','inches','inch','in','mile','mi'
            % SI:  'meters','m','millimeters','mm','kilometers','km'
            %
            % The logic in this class can be found by typing "help
            % conversion" in the command window.
            %
            % This code can quickly be modified to include other units as
            % needed by typing "open convert.m" in the command window.
            
            engToks = {'feet','ft','inches','inch','in','mile','mi'};
            engRat = {1,1,1/12,1/12,1/12,5280,5280};
            engOffset = cell(size(engRat));
            
            siToks = {'kilometers','km','meters','m','centimeters','cm','millimeters','mm'};
            siRat = {1e3,1e3,1,1,1e-2,1e-2,1e-3,1e-3};
            siOffset = cell(size(siRat));
            
            conv = [0.3048,1/0.3048,1]; 
            
            outVal = convert.conversion(value,inUnit,outUnit,engToks,engRat,engOffset,siToks,siRat,siOffset,conv);
        end
        
        function [outVal, outUnit] = temperature(value,inUnit,outUnit)
            %% [outVal, outUnit] = length(value,inUnit,outUnit)
            % This function converts lengths of measurements from the
            % imperial (English) system to the SI system. 
            %
            %
            % ENG -> SI: 1 ft = 0.3048 m = conv
            % SI -> ENG: 1 m = 3.2808 ft = conv
            % SI -> SI or ENG -> ENG:  1 = conv
            % 
            % Supported unit inputs are:
            %
            % ENG: 'feet','ft','inches','inch','in','mile','mi'
            % SI:  'meters','m','millimeters','mm','kilometers','km'
            %
            % The logic in this class can be found by typing "help
            % conversion" in the command window.
            %
            % This code can quickly be modified to include other units as
            % needed by typing "open convert.m" in the command window.
            
            engToks = {'rankine','r','fahrenheit','f',};
            engRat = {1,1,1,1};
            engOffset = {0,0,-459.6699999,-459.6699999};
            
            siToks = {'kelvin','k','celsius','c'};
            siRat = {1,1,1,1};
            siOffset = {0,0,-273.15,-273.15};
            
            conv = [5/9,9/5,1]; 
            
            outVal = convert.conversion(value,inUnit,outUnit,engToks,engRat,engOffset,siToks,siRat,siOffset,conv);
        end
        
        function [outVal, outUnit] = pressure(value,inUnit,outUnit)
            engToks = {'psf','psi','inHg','atm','atmopshere'};
            engRat = {1/2116.8,1/14.7,1/29.9213,1,1};
            engOffset = cell(size(engRat));
            
            siToks = {'Pascal','Pa','kPa','bar','torr','mmHg','atm','atmosphere'};
            siRat = {1/101325,1/101325,1/101.325,1/1.01325,1/760,1/760,1,1};
            siOffset = cell(size(siRat));
            
            conv = [1,1,1]; 
            
            outVal = convert.conversion(value,inUnit,outUnit,engToks,engRat,engOffset,siToks,siRat,siOffset,conv);           
        end
        
        function [outVal, outUnit] = area(value,inUnit,outUnit)
            inUnit = erase(inUnit,'2');
            outUnit = erase(outUnit,'2');

            outVal = convert.length(value,inUnit,outUnit);
            outVal = convert.length(outVal,inUnit,outUnit);
        end
            
        function [outVal, outUnit] = volume(value,inUnit,outUnit)
            inUnit = erase(inUnit,'3');
            outUnit = erase(outUnit,'3');
        
            outVal = convert.length(value,inUnit,outUnit);
            outVal = convert.length(outVal,inUnit,outUnit);
            outVal = convert.length(outVal,inUnit,outUnit);
        end
        
        function [outVal, outUnit] = mass(value,inUnit,outUnit)
            engToks = {'pound','lbm','lb','slug'};
            engRat = {1,1,1,0.031081};
            engOffset = cell(size(engRat));
            
            siToks = {'kg'};
            siRat = {1};
            siOffset = cell(size(siRat));
            
            conv = [0.0453592,1/0.453592,1]; 
            
            outVal = convert.conversion(value,inUnit,outUnit,engToks,engRat,engOffset,siToks,siRat,siOffset,conv); 
        end
        
        function [outVal, outUnit] = conduction(value,inUnit,outUnit)
            engToks = {'BTU/hr-ft-R','BTU/hr-ft-F'};
            engRat = {1,1};
            engOffset = cell(size(engRat));
            
            siToks = {'W/m-K','W/m-C'};
            siRat = {1,1};
            siOffset = cell(size(siRat));
            
            conv = [0.578176,1/0.578176,1]; 
            
            outVal = convert.conversion(value,inUnit,outUnit,engToks,engRat,engOffset,siToks,siRat,siOffset,conv);   
        end
            
        function [outVal, outUnit] = convection(value,inUnit,outUnit)
            engToks = {'BTU/ft2-hr-R','BTU/ft2-hr-F','BTU/ft2-min-R','BTU/ft2-min-F''BTU/ft2-s-R','BTU/ft2-s-F'};
            engRat = {1,1,60,60,3600,3600};
            engOffset = cell(size(engRat));
            
            siToks = {'W/m2-K','W/m2-C'};
            siRat = {1,1};
            siOffset = cell(size(siRat));
            
            conv = [0.176228,1/0.176228,1]; 
            
            outVal = convert.conversion(value,inUnit,outUnit,engToks,engRat,engOffset,siToks,siRat,siOffset,conv);
        end
        
        function [outVal, outUnit] = force(value,inUnit,outUnit)
            engToks = {'lbf'};
            engRat = {1};
            engOffset = cell(size(engRat));
            
            siToks = {'N'};
            siRat = {1};
            siOffset = cell(size(siRat));
            
            conv = [4.44822,1/4.44822,1]; 
            
            outVal = convert.conversion(value,inUnit,outUnit,engToks,engRat,engOffset,siToks,siRat,siOffset,conv);
        end
        
        function [outVal, outUnit] = density(value,inUnit,outUnit)
            engToks = {'lbm/ft3','lbm/in3'};
            engRat = {1, 1/(12^3)};
            engOffset = cell(size(engRat));
            
            siToks = {'kg/m3','g/cm3','g/cc'};
            siRat = {1, 1e-3, 1e-3};
            siOffset = cell(size(siRat));
            
            conv = [16.0185,1/16.0185,1]; 
            
            outVal = convert.conversion(value,inUnit,outUnit,engToks,engRat,engOffset,siToks,siRat,siOffset,conv);
        end 
        
        function [outVal, outUnit] = flow(value,inUnit,outUnit)
            engToks = {'Gallon','Gal','ft3'};
            engRat = {1,1,7.48052;};
            engOffset = cell(size(engRat));
            
            siToks = {'Liter','L','m3'};
            siRat = {1,1,1e3};
            siOffset = cell(size(siRat));
            
            conv = [0.264172, 1/0.264172, 1];
            
            outVal = convert.conversion(value,inUnit,outUnit,engToks,engRat,engOffset,siToks,siRat,siOffset,conv);
        end
        
        function [outVal, outUnit] = heatcapacity(value, inUnit, outUnit)
            engToks = {'BTU/lbm-F'};
            engRat = {1};
            engOffset = cell(size(engRat));
            
            siToks = {'J/kg-K'};
            siRat = {1};
            siOffset = cell(size(siRat));
            
            conv = [1/0.000239, 0.000239, 1];
            
            outVal = convert.conversion(value,inUnit,outUnit,engToks,engRat,engOffset,siToks,siRat,siOffset,conv);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        function outVal = conversion(value,inUnit,outUnit,engToks,engRat,engOffset,siToks,siRat,siOffset,conv)
            %% outVal = conversion(value,inUnit,outUnit,engToks,engRat,siToks,siRat,conv)
            % This function runs all of the actual conversions between
            % units. 
            % [1 inUnit = X stdUnits]
            % [1 outUnit = Y stdUnits]
            %
            % outVal = value * conv * X / Y
            toks = [engToks,siToks]; rat = [engRat,siRat]; offset = [engOffset,siOffset];
            if any(strcmpi(inUnit,toks)) && any(strcmpi(outUnit,toks))
                idxIn = find(strcmpi(inUnit,toks));
                idxOut = find(strcmpi(outUnit,toks));
                if any(strcmpi(inUnit,engToks)) && any(strcmpi(outUnit,siToks))
                    idxConv = 1;
                elseif any(strcmpi(inUnit,siToks)) && any(strcmpi(outUnit,engToks))
                    idxConv = 2;
                else
                    idxConv = 3;
                end
                if isempty(offset)
                    offset{idxIn} = 0;
                    offset{idxOut} = 0;
                else
                    if isempty(offset{idxIn})
                        offset{idxIn} = 0;
                    end
                    if isempty(offset{idxOut})
                        offset{idxOut} = 0;
                    end
                end
                
                outVal = (value - offset{idxIn}) .* conv(idxConv) .* rat{idxIn} ./ rat{idxOut} + offset{idxOut};
          
            else
                error('Incorrect units supplied!')
            end
        end
    end
end

