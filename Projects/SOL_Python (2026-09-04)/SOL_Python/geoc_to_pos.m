function [R, r_geoc] = geoc_to_pos(phi_geoc, lambda, h_geoc, a,b)
% Converts geocentric latitude φ & longitude λ [°] to position vector [x;y;z], 
% optionally at a given altitude h above ellipsoid in the geocentric direction (h_φ)
% Vallado 2013, Eq.(3-10), r_δ & r_K 
% Although not used as often as geod_to_pos, it has certain advantages in some situations.
% → Result is approximate when h_geod (not measured directly along the geocentric radius) is provided instead of h_geoc
%   (but still very good as h_φ ≈ h_ϕ)

% == Sph2Cart(lambda, phi_geoc, a + h_geoc) for b == a

    if nargin < 4
        % Earth = AstroConstants
        a = 6378.137;  % = Earth.r1_km
        b = 6356.75231424518;  % = Earth.r2_km
        if nargin < 3
            h_geoc = 0;
            if nargin < 2
                lambda = 0;
            end
        end
    end

    % in function of φ ; Vallado 2012, Eq.(3-10)
    r_geoc = a ./ sqrt(1 + ((a/b)^2 - 1) * sind(phi_geoc).^2);  % ≠ R_N
    % = a ./ sqrt(1 + e_2 * sind(phi_geoc).^2);  % e_2 = (a/b)^2 - 1;  % using second eccentricity squared e′²
    % = a * sqrt((1 - e2) ./ (1 - e2 * cosd(phi_geoc).^2));  % ← Vallado (3-10)

    % Radius of ellipsoid R(φ) from ellipsoid center to surface (h = 0)
    % r_geoc = GeocentricRadius(phi_geoc, a,b);  
    %        = a ./ sqrt(1 + ((a/b)^2 - 1) .* sind(phi_geoc).^2);
    %        = a ./ sqrt(1 + e2 / (1 - e2) .* sind(phi_geoc).^2);
    %        = a ./ sqrt(1 + e_2 .* sind(phi_geoc).^2);
    %        = a * sqrt(1 ./ (1 + (1/(1 - f)^2 - 1) .* sind(phi_geoc).^2));  

    % r_geoc + h_geoc == Rho

    x = (r_geoc + h_geoc) .* cosd(phi_geoc) .* cosd(lambda);
    y = (r_geoc + h_geoc) .* cosd(phi_geoc) .* sind(lambda);
    z = (r_geoc + h_geoc) .* sind(phi_geoc);

    R = [x; y; z];

end  % function