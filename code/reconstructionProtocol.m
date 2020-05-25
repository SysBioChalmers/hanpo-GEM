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
if ~exist([pwd() '/reconstructionProtocol.m']); error(['Make sure that '...
        'your Current Folder is the one containing the reconstructionProtocol file.']); end
cd ../;  root = [pwd() '/'];
data    = [root 'data/'];
code = [root 'code/'];
cd(code)

%% 3.1 INSTALL RAVEN
%{
Refer to https://github.com/SysBioChalmers/RAVEN/wiki/Installation for
details regarding installation process. Use the pathtool function to add
RAVEN, libSBML, and Gurobi subfolders to MATLAB path.
%}

% CURRENTLY (JULY 2019), SCRIPT REQUIRES RAVEN FROM fix/MILP_tolerance
% BRANCH: https://github.com/SysBioChalmers/RAVEN/tree/fix/MILP_tolerance

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
exportToExcelFormat(model, [root 'scrap/hanpo_step3.3.xlsx']);

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
fid           = fopen([data 'biomass/biomassCuration.csv']);
loadedData    = textscan(fid, '%q %q %q %f','delimiter', ',', 'HeaderLines', 1);
fclose(fid);

BM.name       = loadedData{1};    BM.mets     = loadedData{2};
BM.pseudorxn  = loadedData{3};    BM.coeff    = loadedData{4};

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

% Carbohydrates
indexes = find(contains(BM.pseudorxn, 'carbohydrate'));
equations.mets          = BM.mets(indexes);
equations.stoichCoeffs  = BM.coeff(indexes);
model = changeRxns(model, 'r_4048', equations, 1);

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
clear indexes equations loadedData fid BM biomassRxns ans

%% 3.5 CURATION OF LIPID REACTIONS
% H. polymorpha has unique fatty acid and lipid class compositions. SLIMEr
% explicitly models each lipid moiety, with unique chain distribution, but
% to reduce complexity we will only include a subset of possible chain
% distributions. To do this, files with templates reactions will be
% modified to match the desired chain distributions. First read the file.
fid         = fopen([data '/reconstruction/lipidTemplates.txt']);
loadedData  = textscan(fid, [repmat('%q ', [1, 19]) '%q'], 'delimiter', ...
    '\t', 'HeaderLines',1);
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
cd([code 'lipidMetabolism'])
model = addLipidReactions(template,model,modelSce);

fid         = fopen([data '/reconstruction/lipidTransport.txt']);
loadedData  = textscan(fid, [repmat('%q ', [1, 14]) '%q'], 'delimiter', ...
    '\t', 'HeaderLines', 1);
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
loadedData      = textscan(fid,['%q %q %f' repmat(' %q',[1,13])],...
    'delimiter','\t','HeaderLines',1);
fclose(fid);
template.metName    = loadedData{1};	template.bbID   = loadedData{2};
template.bbMW       = loadedData{3};    template.comps  = loadedData{4};
template.chains = {};
for k = 1:length(loadedData)-4; template.chains(:,k) = loadedData{k+4}; end

model=addSLIMEreactions(template,model,modelSce);
cd(code)

save([root 'scrap/lipids.mat'])
% load([root 'scrap/lipids.mat'])
clear fid loadedData template k toRemove ans
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
modelSce2 = removeReactions(modelSce2,'r_4264',true,true,true);

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
[~,~, addedRxns, model]=fillGaps(model,{modelSce2,modelRhto2},false,true);

% Verify that model can now grow
sol=solveLP(model,1)
printFluxes(model,sol.x)
cd([code 'lipidMetabolism'])
model=scaleLipids(model,'tails');
cd(code)

% Set maximum uptake of carbon source back to 1 mmol/gDCW*hr
model = setParam(model,'lb','r_1808',-1);
model = setParam(model,'lb','r_4041',0);
model = deleteUnusedGenes(model);

% Save workspace
save([root 'scrap/gapfilling.mat'])
% load([root 'scrap/gapfilling.mat'])
clear addedRxns modelSce2 modelRhto2 rxns sol biomassRxns

%% 3.7 SAVE TO GITHUB
% The model is tracked and distributed via a GitHub repository. During the
% reconstruction, one 

% Before saving the model, we will add some extra information.
model.annotation.defaultLB    = -1000; % Default lower bound
model.annotation.defaultUB    = +1000; % Default upper bound
model.annotation.taxonomy     = 'taxonomy/460523';
model.annotation.givenName    = 'Eduard';
model.annotation.familyName   = 'Kerkhoven';
model.annotation.email        = 'eduardk@chalmers.se';
model.annotation.organization = 'Chalmers University of Technology';
model.annotation.note         = 'Draft model accompanying book chapter';
model.id                      = 'hanpo';
model.description             = 'Hansenula polymorpha-GEM';

% Add model information
% Remove sce remnants in subSystems
% As a remnant of the homology based reconstruction, some of the more
% complex grRules have redundancies in subunit configurations. Also remove
% unused metabolites and remove 'sce' prefix from subsystems.
run([code 'curation/cleanupModel']);

newCommit(model);
%% 3.8 SIMULATIONS
% To perform simple FBA, use solveLP
sol = solveLP(model,1);
% Flux distributions can be displayed by printFluxes:
printFluxes(model, sol.x, false);

%% 3.9 MANUAL CURATION
% Modify gene associations of gap-filled reactions. Find reactions
% annotationed with R. toruloides genes.
rhtoRxns = find(contains(model.grRules,'RHTO'));
model.rxnNames(rhtoRxns);
model.grRules(rhtoRxns);

% By BLAST we identify that the gene RHTO_03911 has a homolog in H.
% polymorpha: Hanpo2_15704. We will therefore replace the current grRule
% with the H. polymorpha gene:
model = changeGrRules(model, 'y300009', 'Hanpo2_15704', true);
model = deleteUnusedGenes(model);

% Curate the model so it can support growth on methanol. Reactions required
% for methanol assimilation:
% 1. methanol oxidase:              methanol[p] + oxygen[p] => formaldehyde[p] + hydrogen peroxide[p]
% 2. dihydroxyacetone synthase:     formaldehyde[p] + D-xylulose 5-phosphate[p] => glyceraldehyde 3-phosphate[p] + glycerone[p]
% 3. formate dehydrogenase:         formate[c] + NAD[c] => carbon dioxide[c] + NADH[c]
% 4. peroxisomal catalase:          2 hydrogen peroxide[p] => 2 H2O[p] + oxygen[p]
% 5. dihydroxyacetone kinase:       ATP[c] + glycerone[c] => ADP[c] + dihydroxyacetone phosphate[c] + H+[c]
% Reactions 3-5 are already present in model, manually add reactions 1 & 2:
rxnsToAdd.rxns      = {'MOX', 'DAS'};
rxnsToAdd.equations = {'methanol[p] + oxygen[p] => formaldehyde[p] + hydrogen peroxide[p]',...
'formaldehyde[p] + D-xylulose 5-phosphate[p] => glyceraldehyde 3-phosphate[p] + glycerone[p]'};
rxnsToAdd.rxnNames  = {'methanol oxidase','dihydroxyacetone synthase'};
rxnsToAdd.lb = [0,0];
rxnsToAdd.ub = [1000,1000];
rxnsToAdd.eccodes  = {'1.1.3.13','2.2.1.3'};
rxnsToAdd.grRules  = {'Hanpo2_76277','Hanpo2_95557'};
rxnsToAdd.rxnNotes = {'Methanol metabolism reaction added by manual curation',
'Methanol metabolism reaction added by manual curation'};

model = addRxns(model,rxnsToAdd,3,'',true,true);

% Add methanol (and glucose, was also absent) exchange reactions, and
% diffusion of oxygen and H2O to peroxisome:
model = addRxnsGenesMets(model,modelSce,{'r_4494','r_4391','r_1714','r_1980','r_2098'});
model = setParam(model, 'rev', {'r_4494', 'r_4391'}, 1);
model = setParam(model, 'lb', {'r_4494', 'r_4391'}, -1000);

% Add transport reactions between cytoplasm and peroxisome
model = addTransport(model,'c','p',{'D-xylulose 5-phosphate'...
'methanol', 'glycerone', 'glyceraldehyde 3-phosphate'},true,false,'t_');

% Confirm that growth on methanol is possible
% Ensure glycerol/other carbon source uptake is not allowed
model = setParam(model,'eq',{'r_1808','r_1714'},0);

% Ensure methanol uptake is allowed
model = setParam(model,'lb','r_4494',-1);

% Verify that model can grow
sol=solveLP(model,1)
printFluxes(model,sol.x)

newCommit(model);