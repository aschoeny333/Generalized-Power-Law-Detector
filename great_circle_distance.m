% Calculate the Great Circle Distance Between Receivers
% 
% Author: Alex Schoeny
% 
% Given the latitude and longitude giving the maximum possible distance
% between two receivers, as calculated in possible_range.m, calculate the
% great circle distance between the two receivers (using Haversine
% Formula). Latitude-dependent Earth radius is also calculated.
% 
% Inputs:
%     lat1 - 1 x 1 Double, latitude of reference receiver (in radians)
%     
%     lat2 - 1 x 1 Double, latitude of receiver j (in radians)
%     
%     lon1 - 1 x 1 Double, longitude of reference receiver (in radians)
%     
%     lon2 - 1 x 1 Double, longitude of receiver j (in radians)
%     
% Outputs:
%     dist - 1 x 1 Double, great circle distance between the two receivers
%     (in meters)
%     
%     R - 1 x 1 Double, latitude-dependent Earth radius (in meters)
% 
% Source: https://www.movable-type.co.uk/scripts/gis-faq-5.1.html

function [dist, R, a] = great_circle_distance(lat1, lat2, lon1, lon2)
    % Determine the latitude-dependent Earth radius
    R = (6378 - 21 * sin(lat1)) * 1000;
    
    % Find the differences of latitude and longitude
    dlon = lon2 - lon1;
    dlat = lat2 - lat1;
    
    % Compute distance according to Haversine Formula
    a = sin(dlat/2)^2 + cos(lat1) * cos(lat2) * sin(dlon/2)^2;
    c = 2 * asin(min(1,sqrt(a)));
    dist = R * c;