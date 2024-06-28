% Coordinates
wgs84 = wgs84Ellipsoid('kilometer');

lat = deg2rad(41.3874291);
lon = deg2rad(2.16860934);
alt = 123.456321;

azim = deg2rad(23.36849374);
elev = deg2rad(45.67891234);
range = 12789.0123456;

[x_gs, y_gs, z_gs] = geodetic2ecef(wgs84, lat, lon, alt, "radians");
[x_sat, y_sat, z_sat] = aer2ecef(azim, elev, range, lat, lon, alt, wgs84, "radians");
[e, n, u] = ecef2enu(x_sat, y_sat, z_sat, lat, lon, alt, wgs84, "radians");

fileID = fopen('../data/coordinates.txt', 'w');
fprintf(fileID, 'lat lon alt\n');
fprintf(fileID, '%.15f %.15f %.15f\n', lat, lon, alt);
fprintf(fileID, 'azim elev range\n');
fprintf(fileID, '%.15f %.15f %.15f\n', azim, elev, range);
fprintf(fileID, 'e n u\n');
fprintf(fileID, '%.15f %.15f %.15f\n', e, n, u);
fprintf(fileID, 'x_gs y_gs z_gs\n');
fprintf(fileID, '%.15f %.15f %.15f\n', x_gs, y_gs, z_gs);
fprintf(fileID, 'x_sat y_sat z_sat\n');
fprintf(fileID, '%.15f %.15f %.15f\n', x_sat, y_sat, z_sat);
fclose(fileID);