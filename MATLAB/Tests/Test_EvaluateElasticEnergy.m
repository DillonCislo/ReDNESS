%==========================================================================
% This script is a test of the functionality of 'evaluateElasticEnergy.m'.
% We calculate the elastic energy of a non-Euclidean shell with the
% physical configuration of a flat disk and the target geometry of a
% hemispherical cap.  These results are compared to analytical calculations
%
% By Dillon Cislo 10/4/2019
%==========================================================================

clear; close all; clc;

%--------------------------------------------------------------------------
% Create a triangulation of the flat unit disk
%--------------------------------------------------------------------------
diskTri = diskTriangulationVogel( 1, 2500, 200 );
face = diskTri.ConnectivityList;
v2D = diskTri.Points;

v3D_Disk = [ v2D, zeros( size(v2D,1), 1 ) ];

%--------------------------------------------------------------------------
% Calculate vertex positions on a unit hemispherical cap
%--------------------------------------------------------------------------

r = sqrt( sum( v2D.^2, 2 ) ); % Radial position of each vertex
phi = atan2( v2D(:,2), v2D(:,1) ); % Angular argument of each vertex

thetaCap = pi/2;
v3D_Cap = [ sin( thetaCap .* r ) .* cos(phi) ./ sin(thetaCap), ...
    sin( thetaCap .* r ) .* sin(phi) ./ sin(thetaCap), ...
    ( cos(thetaCap.*r)-cos(thetaCap) ) ./ sin(thetaCap) ];

clear thetaCap r phi

%--------------------------------------------------------------------------
% Calculate the target edge lengths
%--------------------------------------------------------------------------
edgeIDx = diskTri.edges;

tarL = v3D_Cap( edgeIDx(:,2), : ) - v3D_Cap( edgeIDx(:,1), : );
tarL = sqrt( sum( tarL.^2, 2 ) );

%--------------------------------------------------------------------------
% Calculate the target bending angles
%--------------------------------------------------------------------------
resizeCell = @(x) repmat( x, 1, 1+mod(numel(x),2) );
edgeFace = diskTri.edgeAttachments( edgeIDx );
edgeFace = cell2mat( cellfun( resizeCell, edgeFace, ...
    'UniformOutput', false ) );

% Calculate face edges ----------------------------------------------------
e1 = v3D_Cap( face(:,3), : ) - v3D_Cap( face(:,2), : );
e2 = v3D_Cap( face(:,1), : ) - v3D_Cap( face(:,3), : );

% Calculate face unit normals ---------------------------------------------
fN = cross( e1, e2, 2 );
fN =  fN ./ sqrt( sum( fN.^2, 2 ) );

% Calculate bend angles ---------------------------------------------------
tarAngles = dot( fN( edgeFace(:,1), :), fN( edgeFace(:,2), : ), 2 );
tarAngles( tarAngles > 1 ) = 1; tarAngles( tarAngles < -1 ) = -1;
tarAngles = acos( tarAngles );

clear e1 e2 fN edgeFace edgeIDx resizeCell

%--------------------------------------------------------------------------
% Calculate the elastic energy of the configuration
%--------------------------------------------------------------------------

[ E, Es, Eb ] = evaluateElasticEnergy( face, v3D_Disk, tarL, ...
    'TargetAngles', tarAngles );


