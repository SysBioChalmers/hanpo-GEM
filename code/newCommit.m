function newCommit(model)
% newCommit
%   Prepares a new commit of hanpo-GEM in a development branch, for
%   direct submission to GitHub. In contrast to new releases, new commits
%   in development branches do not contain .mat and .xlsx files. This
%   function should be run from the code directory.
%
%   model   RAVEN model structure to be used in commit. If no model
%           structure is provided, the stored XML file will be imported
%           prior to exporting other files
%
%   Usage: newCommit(model)

%Check if in master:
currentBranch = git('rev-parse --abbrev-ref HEAD');
if strcmp(currentBranch,'main')
    error('ERROR: in main branch. For new releases, use newRelease.m')
end

if nargin<1
    %Load model:
    model = importModel('../model/hanpo-GEM.xml');
end

%Save model
exportForGit(model,'hanpo-GEM','../model/',{'txt', 'xml', 'yml'},false,false);
end