function [alignedImage, tform] = imreg2Dr(fixedImage, movingImage)
% imreg2Dr - Perform 2D registration between two grayscale images using intensity-based method
% Syntax:
%   [alignedImage, tform] = imreg2Dr(fixedImage, movingImage)
% Inputs:
%   fixedImage  - Reference image (grayscale)
%   movingImage - Image to be registered/aligned to the fixed image (grayscale)
% Outputs:
%   alignedImage - Registered (aligned) version of the movingImage
%   tform        - Geometric transformation object representing the computed transform
% Description:
%   This function computes a geometric transformation to align movingImage to fixedImage
%   using intensity-based image registration with a translation model by default.
%   It returns the aligned image and the transformation matrix.
% Author:
% Date:


% Define the type of transformation
transformationType = 'translation'; % can be 'translation', 'rigid', 'similarity', 'affine', etc.

% Configure the optimizer and metric for registration
optimizer = registration.optimizer.RegularStepGradientDescent;
metric = registration.metric.MeanSquares;

% Compute the transformation matrix aligning movingImage to fixedImage
tform = imregtform(movingImage, fixedImage, transformationType, optimizer, metric);

% Display the transformation matrix
disp('Transformation Matrix:');
disp(tform.T);

% Apply the geometric transformation to movingImage
alignedImage = imwarp(movingImage, tform, 'OutputView', imref2d(size(fixedImage)));

end
