function [initialDseign detMax] = getMaxsDesigns(initialdesign, dimension,total_size,pointsInDimension, inSigma, inTheta)
% function [initialDseign curMax] = getMaxsDesigns(initialdesign, dimension,total_size,pointsInDimension, inSigma, inTheta)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function initialDseign = getInitialDesignTest()
%
%This is a test function for function getInitialDesign()
%
%PARAMETERS
%dimension - number of design variables
%popSize - number of experiments in the entire design space
%sigma - a prior variance, set to 0.5 in our problem
%theta - a problem dependent parameter, set to 10 in our problem
%
%OUTPUTS
%initialDseign - the initial design obtained using maximum entropy design
%without adaptation
%
%getInitialDesignTest.m     Genzi Li     12/2004
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n\n');
fprintf('You are running the function getInitialDesignTest().\n\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% These parameters can be changed everytime

%dimension = 8; % for test
%initial_size = 600;
%pointsInDimension = 4; %  for test (51*51)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Conduct Some Checks here

if(nargin < 4)
   error('Function must be called with atleast 4 arguments.'); 
end

if(nargin == 4)
   sigma = 0.5;
   theta = 10;
elseif(nargin == 5)
   sigma = inSigma;
   theta = 10;
elseif(nargin == 6)
   sigma = inSigma;
   theta = inTheta;
end

curMax = [];
CovMax = [];


% NEW ONE maxs_augmentOneByOne(currentDesign, dimension, pointsInDimension, designSize, sigma, theta)

% get the initial design using maximum entropy design w/t adaptation
%[initialDseign curMax] = getInitialDesign(dimension, pointsInDimension, initial_size, sigma, theta, initialdesign, Covmax);
[initialDseign, det ]= maxs_augmentOneByOne(initialdesign,  dimension, pointsInDimension, total_size, sigma, theta);

detMax = det;

end % end of function getInitialDesign()
