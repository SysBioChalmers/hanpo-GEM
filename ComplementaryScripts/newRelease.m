function newRelease(bumpType)
% newRelease
%   Prepares a new release of the afmGEM model, for direct submission to
%   GitHub. This function should be run from the ComplementaryScripts
%   directory.
%
%   bumpType    string specifying the type of release, either 'major',
%               'minor' or 'patch'
%
%   Usage: newRelease(bumpType)
%
% Eduard Kerkhoven, 2018-05-15

%Check if in master:
currentBranch = git('rev-parse --abbrev-ref HEAD');
if ~strcmp(currentBranch,'master')
    error('ERROR: not in master')
end

%Bump version number:
oldModel   = load('../ModelFiles/mat/afmGEM.mat');
oldVersion = oldModel.model.description;
oldVersion = oldVersion(strfind(oldVersion,'_v')+2:end);
oldVersion = str2double(strsplit(oldVersion,'.'));
newVersion = oldVersion;
switch bumpType
    case 'major'
        newVersion(1) = newVersion(1) + 1;
        newVersion(2) = 0;
        newVersion(3) = 0;
    case 'minor'
        newVersion(2) = newVersion(2) + 1;
        newVersion(3) = 0;
    case 'patch'
        newVersion(3) = newVersion(3) + 1;
    otherwise
        error('ERROR: invalid input. Use "major", "minor" or "patch"')
end
newVersion = num2str(newVersion,'%d.%d.%d');

%Check if history has been updated:
fid     = fopen('../history.md','r');
history = fscanf(fid,'%s');
fclose(fid);
if ~contains(history,['afmGEMv' newVersion ':'])
    error('ERROR: update history.md first')
end

%Load model:
model = importModel('../ModelFiles/xml/afmGEM.xml');

%Include tag and save model:
model.description = ['Aspergillus fumigatus-GEM_v' newVersion];

%Save model
exportForGit(model,'afmGEM','..',{'mat', 'txt', 'xlsx', 'xml', 'yml'});

%Update version file:
fid = fopen('../version.txt','wt');
fprintf(fid,newVersion);
fclose(fid);
end