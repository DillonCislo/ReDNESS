function [ E, Es, Eb ] = evaluateElasticEnergy( face, vertex, ...
    tarLength, varargin )
%EVALUATEELASTICENERGY This function evaluates the elastic energy of a
%non-Euclidean shell with a given real-space configuration and a specified
%target geometry
%   Detailed explanation goes here

%--------------------------------------------------------------------------
% Set Default Parameter Values
%--------------------------------------------------------------------------

% Material parameters -----------------------------------------------------
h = 0.01;
nu = 0.3;
Y = ( 1 - nu.^2 );

%--------------------------------------------------------------------------
% Input Processing
%--------------------------------------------------------------------------
if(nargin<1), error('Please supply face connectivity list!'); end
if(nargin<2), error('Please supply vertex coordinates!'); end
if(nargin<3), error('Please supply target edge lengths!'); end

% Check the size of the face connectivity list
sizef = size(face);
if ((sizef(2)~=3)||(length(sizef)~=2))
    error('The face list is improperly sized');
end

% Check the size of the vertex list
sizev = size(vertex);
if ((sizev(2)~=3)||(length(sizev)~=2))
    error('The vertex list is improperly sized');
end

% Check if vertex indices exist
if ( max(face(:)) > sizev(1) )
    error('The face list contains an undefined vertex');
elseif ( min(face(:)) < 1 )
    error('The face list contains a vertex index smaller than 1');
end

% Check the size of the edge list
tr = triangulation( face, vertex );
edgeList = tr.edges;
v1 = edgeList(:,1);
v2 = edgeList(:,2);

% Default curvature is plate-like (i.e., flat) ----------------------------
tarTheta = zeros( size( tarLength ) );

% Check the size of the target length list
sizee = size( v1 );
sizetl = size( tarLength );
if ( ~isequal( sizee, sizetl ) )
    tarLength = tarLength';
    sizetl = size(tarLength);
    if ( ~isequal( sizee, sizetl ) )
        error( 'The target length list is improperly sized');
    end
end

for i = 1:length(varargin)
    
    if isa(varargin{i}, 'double')
        continue;
    end
    if isa(varargin{i}, 'logical')
        continue;
    end
    
    % Target geometry parameters ------------------------------------------
    if ~isempty(regexp(varargin{i},'^[Tt]arget[Aa]ngles', 'match'))
        tarTheta = varargin{i+1};
    end
    
    % Material parameters -------------------------------------------------
    if ~isempty(regexp(varargin{i}, '^[Tt]hickness', 'match'))
        h = varargin{i+1};
    end
    if ~isempty(regexp(varargin{i}, '^[Pp]oisson', 'match'))
        nu = varargin{i+1};
    end
    if ~isempty(regexp(varargin{i}, '^[Yy]oung', 'match'))
        Y = varargin{i+1};
    end
    
end

% Process target geometry input -------------------------------------------

% Check the size of the target angle list
sizeta = size( tarTheta );
if ( ~isequal( sizee, sizeta ) )
    tarTheta = tarTheta';
    sizeta = size(tarTheta);
    if ( ~isequal( sizee, sizeta ) )
        error( 'The target angle list is improperly sized');
    end
end

% Update vIDs to match the 0-indexing in C++
% NOTE: This is done internally for the face connectivity list
v1 = v1-1; v2 = v2-1;

%--------------------------------------------------------------------------
% Evaluate Elastic Energy
%--------------------------------------------------------------------------

[ E, Es, Eb ] = evaluate_elastic_energy( face, vertex, tarLength, ...
    tarTheta, int32(v1), int32(v2), h, nu, Y );

end

