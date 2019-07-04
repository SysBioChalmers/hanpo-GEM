%% PROTOCOL FOR THE RECONSTRUCTION OF A Hansenula polymorpha GENOME SCALE MODEL USING THE RAVEN TOOLBOX
%  AUTHORS: FRANCISCO ZORRILLA, EDUARD KERKHOVEN
%
% Complementary to the book chapter XXXX, where more detailed descriptions
% are provided. Section numbering in this script match section numbering
% in the book chapter. 
%
% Before running further code, it can be convenient to define that path
% to some often refered to locations in the model repository. Run these
% lines everytime you start a new MATLAB session where you will be working
% with this script.
clear; clc;
if ~exist([pwd() '\reconstructionProtocol.m']); error(['Make sure that '...
        'your Current Folder is the one containing this Matlab script!']); end
cd ../;  root = [pwd() '/'];
data    = [root 'ComplementaryData/'];
scripts = [root 'ComplementaryScripts/'];
cd(scripts)

%% 3.1 INSTALL RAVEN
%{
Refer to https://github.com/SysBioChalmers/RAVEN/wiki/Installation for
details regarding installation process. Use the pathtool function to add
RAVEN, libSBML, and Gurobi subfolders to MATLAB path.
%}

% Run the following line to ensure that the installation was succesful:
checkInstallation

% This should typically give the following output:
%{
*** THE RAVEN TOOLBOX v. 2.0 ***

Checking if RAVEN is on the Matlab path... PASSED
Checking if it is possible to parse a model in Microsoft Excel format... PASSED
Checking if it is possible to import an SBML model using libSBML... PASSED
Solver found in preferences... gurobi
Checking if it is possible to solve an LP problem using gurobi... PASSED
Checking if it is possible to solve an LP problem using mosek... FAILED
Checking if it is possible to solve an LP problem using cobra... PASSED
Preferred solver... KEPT
Solver saved as preference... gurobi

Checking the uniqueness of RAVEN functions across Matlab path...
No conflicting functions were found
%}

%{
NOTE: Some MATLAB toolboxes may cause conflicitng errors with parsing excel
files. If checkInstallation returns "Checking if it is possible to parse a
model in Microsoft Excel format... FAILED", then please uninstall the
Text Analytics Toolbox. For support for installing and running RAVEN,
please refer to https://gitter.im/SysBioChalmers/RAVEN and 
https://github.com/SysBioChalmers/RAVEN/issues.
%}

%% 3.2 IMPORT TEMPLATE MODELS

% Load S.cerevisiae template GEM using importModel()
% Source: https://github.com/SysBioChalmers/yeast-GEM/releases/tag/v8.3.0
modelSce = importModel([data 'templateModels/yeastGEM.xml']);

% yeast-GEM as obtained from its GitHub repository, is written by the COBRA
% Toolbox. COBRA appends the metabolite compartment to the metabolite ID.
% RAVEN has a separate field (.metComps) for this, so we can remove the
% redundant compartment information from the metabolite IDs.
modelSce.mets = regexprep(modelSce.mets, '\[[a-z]+\]$','');

% Change model ID to 'sce', this is used by getModelFromHomology to match
% each model with its corresponding protein FASTA.
modelSce.id = 'sce';

% If one wants to evaluate what the model contains, it can be easier to
% browse through it as an Excel sheet. During the reconstruction we
% occasionally want to write files that we don't necessarily want to keep
% and track, but just make temporarily. We can write these to a 'scrap'
% folder that will not be tracked on GitHub.
mkdir([root 'scrap'])
exportToExcelFormat(modelSce, [root 'scrap/modelSce.xlsx']);

% Similarly, load R. toruloides template GEM
% Source: https://github.com/SysBioChalmers/rhto-GEM/releases/tab/v1.1.2
modelRhto = importModel([data 'templateModels/rhto.xml'],true);
modelRhto.id = 'rhto';

%{
Note that RAVEN is used to write rhto-GEM in its repository, so there is
no redundant indication of compartments in the metabolite IDs that we
otherwise would have wanted to clear out.
%}

% Ensure that the lower bounds of non-growth associated maintenance and
% biomass formation reactons are set to zero, to make sure we are will not
% force our draft model to produce biomass or regenerate ATP at rates that
% it cannot support.
modelSce  = setParam(modelSce,  'lb', {'r_4041', 'r_4046'}, 0);
modelSce  = setParam(modelSce,  'ub', {'r_4041', 'r_4046'}, 1000);
modelRhto = setParam(modelRhto, 'lb', {'r_4041', 'r_4046'}, 0);
modelRhto = setParam(modelRhto, 'ub', {'r_4041', 'r_4046'}, 1000);

% It can be useful to store intermediate states of the MATLAB environment.
% We will store these in the 'scrap' folder, loading them the next time
% will then allow us to continue from this point.
save([root 'scrap/importModels.mat'])
% load([root 'scrap/importModels.mat'])

%% 3.3 GENERATE MODELS FROM HOMOLOGY
%% 3.3.1 MATCH PROTEIN FASTA IDs
%% 3.3.2 DETERMINE HOMOLOGY by BLAST 

% BLAST the whole-genome protein FASTA of H. polymorpha against the
% S. cerevisiae and R. toruloides protein FASTAs. This can take some
% minutes.
blast = getBlast('hanpo', [data 'genomes/hanpo.faa'], {'sce','rhto'}, ...
    {[data 'genomes/sce.faa'], [data 'genomes/rhto.faa']});
 
% The blast results are then used to generate a new draft model.
model = getModelFromHomology({modelSce, modelRhto}, blast, 'hanpo', ...
    {'sce', 'rhto'}, 1, false, 10^-30, 150, 35);

% Contract model, to merge identical reactions with different grRules
model = contractModel(model);

% Add exchange reactions for media components
mediumComps = {'r_1654', 'r_1672', 'r_1808', 'r_1832', 'r_1861', ...
               'r_1992', 'r_2005', 'r_2060', 'r_2100', 'r_2111'};
model = addRxnsGenesMets(model, modelSce, mediumComps);

% You might want to inspect your first draft model!
exportToExcelFormat(model, [scrap 'hanpo_step3.3.xlsx']);

% Save workspace
save([root 'scrap/homology.mat'])
% load([root 'scrap/homology.mat'])
clear blast mediumComps

%% 3.4 DEFINE BIOMASS COMPOSITION
% H. polymorpha has poly-unsaturated fatty acids, similar to R. toruloides.
% Use the biomass pseudoreactions from rhto-GEM as template to modify.

% Find all reactions with 'pseudreaction' in reactio name in rhto-GEM, and
% add these to the draft model.
biomassRxns = modelRhto.rxns(endsWith(modelRhto.rxnNames, 'pseudoreaction'));
model = addRxnsGenesMets(model, modelRhto, biomassRxns);
% For the lipid curation in step 3.5 and gapfilling in step 3.6 
model = addRxnsGenesMets(model, modelRhto,{'r_4062', 'r_4064', 'r_4046'});
model = setParam(model, 'ub', {'r_4062', 'r_4064', 'r_4046'}, 1000);
model = setParam(model, 'lb', {'r_4062', 'r_4064', 'r_4046'}, 0);

% Load biomass information
fid             = fopen([data 'biomass/biomassCuration.csv']);
loadedData      = textscan(fid, '%q %q %q %f','delimiter', ',', 'HeaderLines', 1);
fclose(fid);

BM.name         = loadedData{1};    BM.mets     = loadedData{2};
BM.pseudorxn    = loadedData{3};    BM.coeff    = loadedData{4};

% Nucleotides (DNA)
% Find out which rows contain the relevant information
indexes = find(contains(BM.pseudorxn, 'DNA'));
% Define new stoichiometries
equations.mets          = BM.mets(indexes);
equations.stoichCoeffs  = BM.coeff(indexes);
% Change reaction
model = changeRxns(model, 'r_4050', equations, 1);

% Ribonucleotides (RNA)
indexes = find(contains(BM.pseudorxn, 'RNA'));
equations.mets          = BM.mets(indexes);
equations.stoichCoeffs  = BM.coeff(indexes);
model = changeRxns(model, 'r_4049', equations, 1);

% Amino acids (protein)
indexes = find(contains(BM.pseudorxn, 'AA'));
equations.mets          = BM.mets(indexes);
equations.stoichCoeffs  = BM.coeff(indexes);
model = changeRxns(model, 'r_4047', equations, 1);

% Lipid backbones
indexes = find(contains(BM.pseudorxn, 'backbone'));
equations.mets          = BM.mets(indexes);
equations.stoichCoeffs  = BM.coeff(indexes);
model = changeRxns(model, 'r_4063', equations, 1);

% Lipid chains
indexes = find(contains(BM.pseudorxn, 'chain'));
equations.mets          = BM.mets(indexes);
equations.stoichCoeffs  = BM.coeff(indexes);
model = changeRxns(model, 'r_4065', equations, 1);

save([root 'scrap/biomass.mat'])
% load([root 'scrap/biomass.mat'])
clear indexes equations loadedData fid BM biomassRxns

%% 3.5 CURATION OF LIPID REACTIONS
% H. polymorpha has unique fatty acid and lipid class compositions. SLIMEr
% explicitly models each lipid moiety, with unique chain distribution, but
% to reduce complexity we will only include a subset of possible chain
% distributions. To do this, files with templates reactions will be
% modified to match the desired chain distributions. First read the file.
fid         = fopen([data '/reconstruction/lipidTemplates.txt']);
loadedData  = textscan(fid, [repmat('%q ', [1, 19]) '%q'], 'delimiter', '\t', 'HeaderLines',1);
fclose(fid);

% Reorganize the content so that it can be used by the addLipidReactions
% function.
template.rxns   = loadedData{1};   template.eqns        = loadedData{2};
template.comps  = loadedData{3};   template.grRules     = loadedData{4};
template.chains = {};
for k = 1:length(loadedData)-4; template.chains(:,k) = loadedData{k+4}; end

% Remove reactions that match a lipid template reaction (ignoring acyl-chains)
toRemove    = regexprep(template.rxns,'CHAIN.*','');
toRemove    = find(startsWith(model.rxnNames,toRemove));
model       = removeReactions(model,toRemove);

% Now use the templates to add the relevant reactions to the model. If a
% reaction already existed in the S. cerevisiae template model, then it
% will use the same reaction identifier.
cd([scripts 'lipidMetabolism'])
model = addLipidReactions(template,model,modelSce);

fid         = fopen([data '/reconstruction/lipidTransport.txt']);
loadedData  = textscan(fid, [repmat('%q ', [1, 14]) '%q'], 'delimiter', '\t', 'HeaderLines', 1);
fclose(fid);
clear template
template.rxns   = loadedData{1};   template.eqns        = loadedData{2};
template.comps  = loadedData{3};   template.chains = {};
for k = 1:length(loadedData)-3; template.chains(:,k) = loadedData{k+3}; end

model = addLipidReactions(template,model,modelSce);

%SLIMER SLIMER SLIMER
clear template
% First remove any SLIME reactions that might exist in the draft model.
model = removeReactions(model,contains(model.rxnNames,'SLIME rxn'));

% Load SLIME template reactions and parse it through the addSLIMEreactions
% function to amend the model.
fid             = fopen([data '/reconstruction/SLIMERtemplates.txt']);
loadedData      = textscan(fid,['%q %q %f' repmat(' %q',[1,13])],'delimiter','\t','HeaderLines',1);
fclose(fid);
template.metName    = loadedData{1};	template.bbID   = loadedData{2};
template.bbMW       = loadedData{3};    template.comps  = loadedData{4};
template.chains = {};
for k = 1:length(loadedData)-4; template.chains(:,k) = loadedData{k+4}; end

model=addSLIMEreactions(template,model,modelSce);
cd(scripts)

save([root 'scrap/lipids.mat'])
% load([root 'scrap/lipids.mat'])
clear fid loadedData template k 
%% 3.6 PERFORM GAP-FILLING

% Use biomass production as obj func for gapfilling
model = setParam(model,'obj','r_4041',1);

% Set biomass production to arbitrary low flux, to force gap-filling to
% produce biomass.
model = setParam(model,'lb','r_4041',0.01);

% Set glycerol uptake at a higher value, to make sure that fluxes are high
% enough during gap-filling that they won't be ignored due to the tolerance
% of the MILP solver
model = setParam(model,'lb','r_1808',-10);

% From the Sce and Rhto models, remove all exchange reactions (the
% necessary ones we already added, don't want to add new ones)
modelSce2 = removeReactions(modelSce,getExchangeRxns(modelSce,'both'),true,true,true);
biomassRxns = modelSce2.rxns(endsWith(modelSce2.rxnNames,'pseudoreaction'));
modelSce2 = removeReactions(modelSce2,biomassRxns,true,true,true);

modelRhto2 = removeReactions(modelRhto,getExchangeRxns(modelRhto,'both'),true,true,true);
% Also, find which reactions in Rhto are already present in Sce (according
% to reaction ID, as Rhto model is based on Sce). No need to have duplicate
% reactions as suggestions during gap-filling, this reduces the solution
% space and helps to solve the MILP part of gap-filling.
rxns = intersect(modelSce2.rxns,modelRhto2.rxns);
modelRhto2 = removeReactions(modelRhto2,rxns);
biomassRxns = modelRhto2.rxns(endsWith(modelRhto2.rxnNames,'pseudoreaction'));
modelRhto2 = removeReactions(modelRhto2,biomassRxns,true,true,true);

% Run fillGaps function
[~,~, addedRxns, ~, ~]=fillGaps(model,{modelSce2,modelRhto2},false,true);

SceRxns = addedRxns(ismember(addedRxns,modelSce.rxns));
allRxns = getAllRxnsFromGenes(modelSce,SceRxns);
model   = addRxnsGenesMets(model, modelSce, allRxns, true);
RhtoRxns = addedRxns(ismember(addedRxns,modelRhto2.rxns))
allRxns  = getAllRxnsFromGenes(modelRhto,RhtoRxns);
model    = addRxnsGenesMets(model, modelRhto, allRxns, true);

% Verify that model can now grow
sol=solveLP(model,1)
printFluxes(model,sol.x)
cd([scripts 'lipidMetabolism'])
model=scaleLipids(model,'tails');

% Set maximum uptake of carbon source back to 1 mmol/gDCW*hr
model = setParam(model,'lb','r_1808',-1);

% Save workspace
save([root 'scrap/gapfilling.mat'])
% load([root 'scrap/gapfilling.mat'])
clear addedRxns modelSce2 modelRhto2 rxns sol allRxns RhtoRxns SceRxns

%% 3.7 MANUAL CURATION

% Replace gap-filled genes with orthologs

% Curate some other reactions

% Save workspace
save([root 'scrap/curation.mat'])
% load([root 'scrap/curation.mat'])

%% 3.6 SAVE TO GITHUB


% Add model information
% Remove sce remnants in subSystems

% https://github.com/SysBioChalmers/Hansenula_polymorpha-GEM

% Create github repository similar to one linked above

% Upload first homologyModel draft

%% 3.7 SIMULATIONS


