function runDemoDstorm(srcDir)
%RUNDEMODSTORM Run the 3D-dSTORM septum-width module on a folder of profile CSVs.
%
%   runDemoDstorm()          uses demo/demo_output/dSTORM_profiles
%   runDemoDstorm(srcDir)    uses the given folder of intensity-profile CSVs
%
%   The wrapper exists so that the module script dSTORMwidthcaclu.m (which
%   starts with "clearvars" and "clc") runs in its own workspace and does not
%   clear the variables of the calling demo script.

if nargin < 1 || isempty(srcDir)
    srcDir = fullfile(fileparts(mfilename('fullpath')), 'demo_output', 'dSTORM_profiles');
end

repoDir = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repoDir, '3D-dSTORM'));

sourceFolder = srcDir;      % consumed by dSTORMwidthcaclu.m
run(fullfile(repoDir, '3D-dSTORM', 'dSTORMwidthcaclu.m'));

end
