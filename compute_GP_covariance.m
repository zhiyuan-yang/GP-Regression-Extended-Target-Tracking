function [ covMatrix ] = compute_GP_covariance(argArray1, argArray2, sigmaf,sigmar,l)
% This function computes the GP covariance according to the specified parameters.
% The type of the kernel function is the squared exponential. 
% Output:
%           covMatrix:      The covariance matrix computed by the GP kernel.
% Inputs:
%           argArray1:      First argument of the covariance function. (Column vector)
%           argArray2:      Second argument of the covariance function. (Column vector)
%           paramGP:        Parameters of the specified GP


% Author:   Murat Kumru <kumru@metu.edu.tr>


% Extract the parameters
stdPrior = sigmaf;
stdRadius  = sigmar;
lengthScale = l;

diffMatrix = compute_diffence_matrix(argArray1, argArray2);

% Covariance periodic with 2*pi
covMatrix = stdPrior^2 * exp(-2 * sin(diffMatrix/2).^2 / lengthScale^2) + stdRadius^2;

end


function [ diffMatrix ] = compute_diffence_matrix(argArray1, argArray2)

len1 = length(argArray1);
len2 = length(argArray2);

argGrid1 = repmat(argArray1, [1, len2]);
argGrid2 = repmat(transpose(argArray2), [len1, 1]);

diffMatrix = argGrid1-argGrid2;
end
