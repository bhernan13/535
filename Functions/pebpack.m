function [npeb_layer, porosity, pebfill, vol_total] = pebpack(dpeb, dpipe)

%% Pebble packing calculation
% Created: February 16, 2019
% Code Author: Madie Melcer
% Code to calculate packing of pebbles for one layer of pebble bed given
% pebble diameter and bed inner diameter. Assumes that layers of pebbles do
% not overlap and that hexagonal packing of pebbles is used (area fill
% factor of pi/(12^0.5).
% Variables
% dpeb = pebble diameter
% dpipe = diameter pipe
% npeb_layer = pebbles per layer - circles of diameter dpeb fit into circle
% of diameter dpipe
% vol_peb = total volume of the pebbles
% vol_total = volume of the cylinder of diameter dpipe and height equal to
% dpeb
% vol_void = volume of the unfilled space 
% porosity = ratio of void volume to total volume
% pebfil = ratio of pebble volume to total volume
npeb_layer = floor((dpipe./dpeb).^2*pi./sqrt(12));

vol_peb = npeb_layer.*( ((4/3).*pi).*((dpeb./2).^3) ); % volume of pebbles in one layer
vol_total = (pi/4).*(dpipe.^2).*dpeb; % volume for slice of cylinder with height 
    % equal to dpeb
vol_void = vol_total-vol_peb;
porosity = vol_void./vol_total;
pebfill = vol_peb./vol_total;