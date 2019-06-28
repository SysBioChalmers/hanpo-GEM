%% PROTOCOL FOR THE RECONSTRUCTION OF A Hansenula polymorpha GENOME SCALE METABOLIC MODEL USING THE RAVEN TOOLBOX
%  AUTHORS: FRANCISCO ZORRILLA, EDUARD KERKHOVEN

%3. Method
%   3.1 INSTALL RAVEN
%   3.2 IMPORT TEMPLATE MODELS
%       3.2.1 SCE
%       3.2.2 RHTO
%   3.3 GENERATE TARGET MODELS FROM HOMOLOGY
%       3.3.1 SCE
%       3.3.2 RHTO
%   3.4 MERGE MODELS
%   3.5 GENERATE BOF
%   3.6 GAP-FILLING
%   3.7 SAVE TO GITHUB
%   3.8 SIMULATIONS

%% 3.1 INSTALL RAVEN

%{

Refer to https://github.com/SysBioChalmers/RAVEN/wiki/Installation for
details regarding installation process. Use the pathtool function to add
RAVEN, COBRA, and Gurobi subfolders to MATLAB path.

%}

% Run the following line to ensure that the installation was succesful:
checkInstallation

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
Text Analytics Toolbox.

%}

%% 3.2 IMPORT TEMPLATE MODELS

cd C:\Users\zorrilla\Desktop\GEM\bookChapter

%% 3.2.1 IMPORT S.cerevisiae MODEL

% Load S.cerevisiae template GEM using importModel()
% Source: https://github.com/SysBioChalmers/yeast-GEM/releases/tag/v8.3.0
modelSce=importModel('yeastGEMv8.3.0.xml');

% Remove compartment information from metabolite ids (redundant).
modelSce.mets=regexprep(modelSce.mets,'\[[a-z]+\]$','');

% Change model ID to "sce" for clarity
modelSce.id = 'sce';

% Change reversibility fields to match boundaries. Prevents problems with
% MENECO gapfilling.
for i=1:length(modelSce.rxns)
	if modelSce.lb(i)==0 && not(modelSce.ub(i)==0);
        modelSce.rev(i)=0;
    elseif modelSce.lb(i)<=0 && modelSce.ub(i)>=0;
        modelSce.rev(i)=1;
	end
end

% Set the objective function of the model to growth using setParam() and 
% confirm that the model is functional by using solveLP(). Reaction 2111
% correponds to growth in this particular model, this can be double 
% checked in the excel file.
modelSce=setParam(modelSce,'obj','r_2111',1);

% Perform Flux Balance Analysis using solveLP()
sol=solveLP(modelSce,1)

% Print exchange fluxes
printFluxes(modelSce,sol.x)

% Export model as excel format for user-friendly inspection
exportToExcelFormat(modelSce,'modelSce.xlsx');

% Export model as .xml format to be used for gapfilling
exportModel(modelSce, 'modelSce.xml');

clear i sol

%% 3.2.2 IMPORT R.toruloides MODEL

% Load R.toruloides template GEM using importModel()
% Source: https://github.com/SysBioChalmers/rhto-GEM/blob/master/ModelFiles/xml/rhto.xml

modelRhto=importModel('rhto.xml',true);

for i=1:length(modelRhto.rxns)
	if modelRhto.lb(i)==0 && not(modelRhto.ub(i)==0);
        modelRhto.rev(i)=0;
    elseif modelRhto.lb(i)<=0 && modelRhto.ub(i)>=0;
        modelRhto.rev(i)=1;
	end
end

% Set the objective function to growth
modelRhto=setParam(modelRhto,'obj','r_2111',1);

% Perform Flux Balance Analysis using solveLP()
sol=solveLP(modelRhto,1)

% Print exchange fluxes
printFluxes(modelRhto,sol.x)

% Export model to excel for inspection
exportToExcelFormat(modelRhto,'modelRhto.xlsx');

clear sol i

% Save workspace
save('C:\Users\zorrilla\Desktop\GEM\bookChapter\IMPORTMODELS.mat')
% load('C:\Users\zorrilla\Desktop\GEM\bookChapter\IMPORTMODELS.mat')

%% 3. CONSTRUCT DRAFT MODEL FROM HOMOLOGY

%{

To construct a draft model from homology we can use the protein FASTA 
files obtained in the PREPARATION section. Make sure gene identifiers match
as mentioned previously!

First, we perform a biderectional protein BLAST search using the getBlast()
function, to determine homologous proteins/enzyme catalyzed reactions in 
H.polymorpha and S.cerevisiae. This step results in a BLAST strucutre.


%}

%% 3.3.1 GENERATE TARGET MODEL FROM HOMOLOGY USING S.cerevisiae MODEL

%hpoSceBlast = getBlast('hpo','hpo.faa',{'sce'},{'sce.faa'});

% Save the BLAST structure to avoid unnecessarily re-running getBlast() in
% a future session
%save('hpoSceBlast.mat','hpoSceBlast');

% Next time use load() to access the BLAST structure
load('hpoSceBlast.mat');

% Obtain model for H. polymorpha using BLAST structure and the expanded
% S. cerevisiae model using getModelFromHomology().
modelHpoHomology1=getModelFromHomology({modelSce},hpoSceBlast,'hpo',{},1,false,10^-20,100,35,true);

% Add all exchange reactions. These were not gene annotated, and therefore 
% not added by the getModelFromHomology() function. We might not require 
% all exchange rxns, but we can easily remove unconnected ones later.
exRxns=getExchangeRxns(modelSce);                                                
modelHpoHomology1=addRxnsGenesMets(modelHpoHomology1,modelSce,exRxns,false,'Modeling reaction',1);

% Same as exchange reactions, add all non-gene annotated transport reactions
noGeneIdx=find(cellfun(@isempty,modelSce.grRules)); % Find rxns with no genes
rxnIdx=regexp(modelSce.rxnNames,'(transport)|(diffusion)|(pseudoreaction)'); % Find rxnNames containing 'transport' or 'diffusion' or 'pseudoreaction'
rxnIdx=find(~cellfun('isempty',rxnIdx)); % Find their indices
rxnIdx=intersect(rxnIdx,noGeneIdx); % Keep the ones without gene notation
exRxns=modelSce.rxns(rxnIdx); % Obtain reaction IDs
modelHpoHomology1=addRxnsGenesMets(modelHpoHomology1,modelSce,exRxns,false,...
    'Modeling reaction required for intercellular transport, gene unknown',1); % Finally add reactions and metabolites

% Export model to excel for inspection
exportToExcelFormat(modelHpoHomology1,'modelHpoHomology1.xlsx'); 

% Clean up workspace
clear rxnIdx noGeneIdx exRxns modelSceExpanded hpoSceBlast sol

%% 3.3.2 GENERATE TARGET MODEL FROM HOMOLOGY USING R.toruloides MODEL

%hpoRhtoBlast = getBlast('hpo','hpo.faa',{'rhto'},{'rhto_np11.faa'});

% Save the BLAST structure to avoid unnecessarily re-running getBlast() in
% a future session
%save('hpoRhtoBlast.mat','hpoRhtoBlast');

% Next time use load() to access the BLAST structure
load('hpoRhtoBlast.mat');

% Obtain model for H. polymorpha using BLAST structure and the expanded
% R. toruloides model using getModelFromHomology().
modelHpoHomology2=getModelFromHomology({modelRhto},hpoRhtoBlast,'hpo',{},1,false,10^-20,100,35,true);

% Export model to excel for inspection
exportToExcelFormat(modelHpoHomology2,'modelHpoHomology2.xlsx'); 

clear exRxns hpoRhtoBlast modelRhtoExpanded noGeneIdx rxnIdx sol

% Save workspace
save('C:\Users\zorrilla\Desktop\GEM\bookChapter\HOMOLOGY.mat')
% load('C:\Users\zorrilla\Desktop\GEM\bookChapter\HOMOLOGY.mat')

%% 3.4 MERGE HOMOLOGY-DERIVED MODELS

% Use mergeModels function to merge the two homology-derived drafts
modelHpoMerged = mergeModels({modelHpoHomology1,modelHpoHomology2})

% Contract obtained model
modelHpoMerged = contractModel(modelHpoMerged)

% Change ID from MERGED to hpo
modelHpoMerged.id = 'hpo';

% Add model description
modelHpoMerged.description = 'H. polymorpha GEM genearated by merging 2 homology based models from S.ce and Rh.to';

% Set objective function to growth (biomass equation taken from S.ce
% template)
modelHpoMerged=setParam(modelHpoMerged,'obj','r_2111',1);

% Verify that model cannot grow, most likely due to gaps
sol=solveLP(modelHpoMerged,1)
printFluxes(modelHpoMerged,sol.x)

% Export model as a .xml file
exportModel(modelHpoMerged, 'modelHpoMerged.xml');

% Export model to excel for inspection
exportToExcelFormat(modelHpoMerged,'modelHpoMerged.xlsx'); 

clear sol

% Save workspace
save('C:\Users\zorrilla\Desktop\GEM\bookChapter\MERGEMODELS.mat')
% load('C:\Users\zorrilla\Desktop\GEM\bookChapter\MERGEMODELS.mat')

%% 3.5 GENERATE BIOMASS OBJECTIVE FUNCTION COEFFICIENTS

% Rename dNTPs metabolite IDs to bigg ID for BOFdat:


%BOFdat output:

%>>> dna_coefficients

%{<Metabolite datp_c at 0x7f2b6904fd90>: -0.00855716169972182, 
%<Metabolite dctp_c at 0x7f2b6904fe90>: -0.008628417496747516, 
%<Metabolite ppi_c at 0x7f2b690aa250>: 0.03354733319457753, 
%<Metabolite dgtp_c at 0x7f2b69064250>: -0.007579394372551919, 
%<Metabolite dttp_c at 0x7f2b69064910>: -0.00878235962555628}

rxnIndex= getIndexes(modelHpoMerged,{'r_4050'},'rxns')
modelHpoMerged.S(:,index)
modelHpoMerged.S(380,index)= -0.00855716169972182;
modelHpoMerged.S(385,index)= -0.008628417496747516;
modelHpoMerged.S(397,index)= -0.007579394372551919;
modelHpoMerged.S(419,index)= -0.00878235962555628;

clear rxnIndex

exportToExcelFormat(modelHpoMerged,'modelHpoMerged.xlsx'); 

% Save workspace
save('C:\Users\zorrilla\Desktop\GEM\bookChapter\BIOMASS.mat')
% load('C:\Users\zorrilla\Desktop\GEM\bookChapter\MERGEMODELS.mat')


%% 3.6 PERFORM GAP-FILLING

% Find targets: any substrate for the pseudoreactions, use the following
% text to reconstruct menecoTargets.sbml.
rxnIdx=find(contains(modelHpoMerged.rxnNames,'pseudoreaction'));
targets=find(any(modelHpoMerged.S(:,rxnIdx)<0,2));
[modelHpoMerged.mets(targets), modelHpoMerged.metNames(targets)]
targetSBML=strcat('<species id="M_',modelHpoMerged.mets(targets),...
     '" name="',modelHpoMerged.metNames(targets),'"/>');

% Run MENECO outside of MATLAB
 
% Import identified reactions
fid     = fopen ('C:\Users\zorrilla\Desktop\GEM\bookChapter\rxns.txt');
menecoRxns = textscan(fid,'%s'); fclose(fid);
menecoRxns = menecoRxns{1};

% Add reactions
menecoRxns=getAllRxnsFromGenes(modelSce,menecoRxns);
modelHpoGF=addRxnsGenesMets(modelHpoMerged,modelSce,menecoRxns,true,...
'Identified by MENECO to produce biomass components',1); % Add reactions and metabolites

% Model still cannot grow because mannan is not able to be produced, use
% RAVEN fillGaps function

% Force model to push flux through growth reaction by setting the lower 
% bound to an arbitrary positive number
modelHpoGF = setParam(modelHpoGF,'lb','r_4041',0.01);
% Run fillGaps() function to fill gaps in the model based on growth
% constraint
[newConnected, cannotConnect, addedRxns, modelHpoGF, exitFlag]=fillGaps(modelHpoGF,modelSce,false,true);
% Reset lower bound of growth reaction to zero in gap filled model to test
% if model can now grow
modelHpoGF = setParam(modelHpoGF,'lb','r_4041',0);

% Ensure objective function set to growth
modelHpoGF=setParam(modelHpoGF,'obj','r_4041',1);
% Solve for fluxes
sol=solveLP(modelHpoGF,1);
% Check fluxes
printFluxes(modelHpoGF,sol.x); % Growth achieved

clear addedRxns ans cannotConnect exitFlag fid menecoRxns modelHpoHomology1 modelHpoHomology2 modelHpoMerged newConnected sol targetSBML targets rxnIndex

% Export model to excel for inspection
exportToExcelFormat(modelHpoGF,'modelHpoGF.xlsx'); 

% Export model as a .xml file
exportModel(modelHpoGF, 'modelHpoGF.xml');

% Save workspace
save('C:\Users\zorrilla\Desktop\GEM\bookChapter\GAPFILLING.mat')
% load('C:\Users\zorrilla\Desktop\GEM\bookChapter\GAPFILLING.mat')

%% 3.7 SAVE TO GITHUB

% https://github.com/SysBioChalmers/Hansenula_polymorpha-GEM

% Create github repository similar to one linked above

% Upload first homologyModel draft

%% 5. CONSTRUCT DRAFT MODELS FROM KEGG 

%{

In this section we make use of the KEGG-based resconstruction approaches
offered by the RAVEN toolbox. The first approach is based on hidden markov 
models (HHMs), while the second approach is based on annotation.

%}

% Create model using KEGG & HMM
modelAfmKeggHMM = getKEGGModelForOrganism('afm','afm_prot.faa','euk100_kegg82','output',false,false,false,10^-30,0.8,0.3,-1)
   
exportToExcelFormat(modeAfmKeggHMM,'modelAfmKeggHMM.xlsx') % Inspect model

% Create model using KEGG & annotation
modelAfmKeggAnnot =getKEGGModelForOrganism('afm','','','',false,false,false)   %check if this worked corrctly

exportToExcelFormat(modeAfmKeggAnnot,'modelAfmKeggAnnot.xlsx') % Inspect model

%% 6. CONSTRUCT DRAFT MODEL FROM METACYC
%
modelAfmMC=getMetaCycModelForOrganism('afm','afm_prot.faa',false,false,false,100,55) 

exportToExcelFormat(modelAfmMC,'modelAfmMC.xlsx');

%% 7. COMBINE MODELS 
% 

%{

This is perhaps the trickiest part of the reconstruction procedure, as
merging models can be extremely tedious and time consuming if done in a 
manual fashion. On the other hand, taking an automated approach may result
in many errors when matching metabolites and/or reactions. One approach 
which has been developed to tackle this challange is the modelBorgifier 
Toolbox, which is part of the COBRA toolbox.

%}

% Use mergeModels() function to combine the two KEGG-based approach models
modelAfmKegg= mergeModels({modelAfmKeggHMM,modelAfmKeggAnnot}) %these models have missing names for reactions

%{

modelAfmKegg = 

  struct with fields:

             id: 'MERGED'
    description: ''
           rxns: {3246×1 cell}
       rxnNames: {3246×1 cell}
        eccodes: {3246×1 cell}
     subSystems: {3246×1 cell}
     rxnMiriams: {3246×1 cell}
              S: [1886×3246 double]
           mets: {1886×1 cell}
            rev: [3246×1 double]
             ub: [3246×1 double]
             lb: [3246×1 double]
              c: [3246×1 double]
              b: [1886×1 double]
          genes: {2774×1 cell}
     rxnGeneMat: [3246×2774 double]
       metNames: {1886×1 cell}
    metFormulas: {1886×1 cell}
         inchis: {1886×1 cell}
     metMiriams: {1886×1 cell}
          comps: {'s'}
      compNames: {'System'}
    compOutside: {''}
       metComps: [1886×1 double]
        grRules: {3246×1 cell}
        rxnFrom: {3246×1 cell}
        metFrom: {1886×1 cell}
       geneFrom: {2774×1 cell}
    geneMiriams: {2774×1 cell}

%}

exportToExcelFormat(modelAfmKegg,'modelAfmKegg.xlsx');

%use combineMetaCycKEGGModels() instead, NOTE these dont have compartmens

modelAfmMCKegg = combineMetaCycKEGGModels(modelAfmMC55,modelAfmKegg)

exportToExcelFormat(modelAfmMCKegg,'modelAfmMCKegg.xlsx');

save(modelAfmMCKegg,'modelAfmMCKegg.mat')
save('modelAfmHomology.mat')

%% Localization

%{

This aspect of GEM reconstruction can also be tricky. Since most draft
models obtained up to this point lack any type of localization information
(e.g modelAfmKegg,modelAfmMC, modelAfmMCKegg, etc), it may be best to merge
models first and then predict localization on the final model. The first
step is to break our target GEM organism's FASTA file into chunks of 1000
protein sequences. These chunks can then be run through the subcellular 
localization prediction tool CELLO (http://cello.life.nctu.edu.tw/). The
output of this tool can then be reassembled into one file containing the
localization information for the entire FASTA file. This is then used as an
input for the parseScores() function, which yields a Gene Scoring
Structure (GSS) as the output. This GSS is then subsetted to only include
localization information about the compartments of interest

%}

    % remove localization in homology mergeCompartments() (keep unconstrained), merge models,gap fill ,predict localization
    %use homology/pch modelsz for gap filling 
    
    %do not want all of the compartments from cello, use copyToComps() 
    %function to move rxns arounnd/then delete unwanted comps (or
    %simply subset the obtained GSS to only include compartments of
    %interest). keep cytosol,mitochondria, peroxisome, extracellular; 
    %everything else is cytosolic, what about boundary?
    
    %DO:
    %merge compartments in AfmHomology/Pch, fillGaps, predict localization
    %(compare to original AfmHomology to see if predict localization is working), 
    % then rerun this starting with AfmHomologyKeggMC
    
    %predict localization on KEGG, MC models, if works merge with Homology

GSS = parseScores('afm_prot_deeploc.txt','deeploc'); 
GSSsub = GSS.scores(:,[3 4 5 11]); %only care about extracellular, cytoplasmic, mitochondrial, and peroxisomal
GSSsubset = GSS;
GSSsubset.scores = GSSsub;
GSSsubset.compartments = GSS.compartments([3 4 5 11]);

GSSsubset.compartments = transpose(GSSsubset.compartments) %need to do this for predictLocalization to work

%[modelAfmHomologyMergeComps, deletedRxnsMergeComps, duplicateRxnsMergeComps]=mergeCompartments(modelAfmHomology,true);

%fill gaps before localization, merge comps of Pch model before fill gaps

%[modelPchMerged, deletedRxns, duplicateRxns]=mergeCompartments(modelPch,true);

[newConnected, cannotConnect, addedRxns, newModel, exitFlag]=fillGaps(modelAfmHomologyMergeComps,{modelPchMerged})

[modelAfmHomologyDeeploc, geneLocalization, transportStruct, score, removedRxns] = predictLocalization(modelAfmHomology,GSSsubset,'Cytoplasm')

modelAfmHomologyDeeploc = rmfield(modelAfmHomologyDeeploc,'geneComps');
exportToExcelFormat(modelAfmHomologyDeeploc,'modelAfmHomologyDeeploc.xlsx'); 

% try contracting model, gap
% fill using unmergedPch (careful with naming)

% try same thing but with the modelAfmHomologyPreMerged

%[modelAfmHomologyPreMergedCello, geneLocalization, transportStruct, score, removedRxns] = predictLocalization(modelAfmHomologyPreMerged ,GSSsubset,'cytop')

% localization is kinda jank, try using wolf localization on afm first,
% then on pch to get sense of how well it performs
% check (1.3.3.6 E.C., beta ox in peroxisome and ppp in cytosol)
%otherwise try CELLO on pch 

%% 8. BIOMASS EQUATION checkout BOF dat

%{
%look at pch excel model, there is tab for biomass sc

No hits for biomass composition of A.fumigatus, use A.niger biomass
composition from
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2290933/'



bofDAT (python) protocol:

    1. Obtain the macromolecular composition of the organism

    2. Obtain OMICs experimental data

    3. Generate BOFsc

    4. Generate NGAM and GAM

    5. Update BOF (BOFdat!)


%}


%% 9. GAP FILLING 
%

%% 10. SIMULATIONS
%

%% 11. CURATION
%

%% Final clean up 

% Delete unused genes
% Model information
modelCal2=deleteUnusedGenes(modelCal2);

modelCal2.id='C.alb';
modelCal2.description='Genome-scale model of Candida albicans';
modelCal2.annotation.familyName='Zorrilla';
modelCal2.annotation.givenName='Francisco';
modelCal2.annotation.email='zorrilla@chalmers.se';
modelCal2.annotation.organization='Chalmers University of Technology';
%modelCal2.annotation.taxonomy='4953';
modelCal2.annotation.notes='This was an amateur attempt at creating a GEM, please do not cite';

