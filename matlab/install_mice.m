matlabrc; clc; close all;

% Determine the correct operating system:
if ismac
    error('MAC IS NOT YET SUPPORTED')
elseif isunix
    error('UNIX IS NOT YET SUPPORTED')
elseif ispc
    % Download JPL's MICE:
    url = 'https://naif.jpl.nasa.gov/pub/naif/toolkit//MATLAB/PC_Windows_VisualC_MATLAB9.x_64bit/packages/mice.zip';
    filename = 'mice.zip';
    websave(filename,url);
    
    % Unzip MICE:
    unzip(filename,'.')
    
    % Cleanup:
    delete *.zip
    delete *.html
else
    disp('Platform not supported')
end