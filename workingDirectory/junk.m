
hpoSceBlast = getBlast('hpo','hpo.faa',{'sce'},{'sce.faa'});

% Save the BLAST structure to avoid unnecessarily re-running getBlast() in
% a future session
save('hpoSceBlast.mat','hpoSceBlast');

% Next time use load() to access the BLAST structure
load('hpoSceBlast.mat');

% Obtain draft model for H. polymorpha using BLAST structure and S. cerevisiae 
% model template using getModelFromHomology()
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
clear rxnIdx noGeneIdx exRxns hpoSceBlast

%% 3.3.3 GENERATE TARGET MODEL FROM HOMOLOGY USING R.toruloides MODEL

hpoRhtoBlast = getBlast('hpo','hpo.faa',{'rhto'},{'rhto_np11.faa'});

% Save the BLAST structure to avoid unnecessarily re-running getBlast() in
% a future session
save('hpoRhtoBlast.mat','hpoRhtoBlast');

% Next time use load() to access the BLAST structure
load('hpoRhtoBlast.mat');

% Obtain model for H. polymorpha using BLAST structure and the expanded
% R. toruloides model using getModelFromHomology().
modelHpoHomology2=getModelFromHomology({modelRhto},hpoRhtoBlast,'hpo',{},1,false,10^-20,100,35,true);

% Export model to excel for inspection
exportToExcelFormat(modelHpoHomology2,'modelHpoHomology2.xlsx'); 

clear exRxns hpoRhtoBlast noGeneIdx rxnIdx

% Save workspace
save('HOMOLOGY.mat')
% load('HOMOLOGY.mat')

%% 3.4 MERGE HOMOLOGY-DERIVED MODELS

% Use mergeModels function to merge the two homology-derived drafts
modelHpoMerged = mergeModels({modelHpoHomology1,modelHpoHomology2})

% Contract obtained model
modelHpoMerged = contractModel(modelHpoMerged)

% Change ID from MERGED to hpo
modelHpoMerged.id = 'hpo';

% Add model description
modelHpoMerged.description = 'H. polymorpha GEM genearated by merging two homology based models from S.ce and Rh.to';

% Export model as a .xml file
exportModel(modelHpoMerged, 'modelHpoMerged.xml');

% Export model to excel for inspection
exportToExcelFormat(modelHpoMerged,'modelHpoMerged.xlsx'); 

% Save workspace
save('MERGEMODELS.mat')
% load('MERGEMODELS.mat')

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

