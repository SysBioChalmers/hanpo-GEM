%% PROTOCOL FOR THE RECONSTRUCTION OF A Hansenula polymorpha GENOME SCALE METABOLIC MODEL USING THE RAVEN TOOLBOX
%  AUTHORS: FRANCISCO ZORRILLA, EDUARD KERKHOVEN

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

% Load S.cerevisiae template GEM using importModel()
% Source: https://github.com/SysBioChalmers/yeast-GEM/releases/tag/v8.3.0
modelSce=importModel('yeastGEM.xml');

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

% Export model as excel format for user-friendly inspection
exportToExcelFormat(modelSce,'modelSce.xlsx');

clear i 

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

% Set lb of growth and NGAM to zero, to make sure we're not trying to force
% the model to produce biomass or NGAM at a rate it can't support.
modelSce = setParam(modelSce,'lb',{'r_4041','r_4046'},0);
modelRhto = setParam(modelRhto,'lb',{'r_4041','r_4046'},0);

% Export model to excel for inspection
exportToExcelFormat(modelRhto,'modelRhto.xlsx');

clear i

% Save workspace
save('IMPORTMODELS.mat')
% load('IMPORTMODELS.mat')

%% 3.3 CONSTRUCT DRAFT MODEL FROM HOMOLOGY

%% 3.3.1 MATCH PROTEIN FASTA IDs

%% 3.3.2 GENERATE TARGET MODEL FROM HOMOLOGY 

% Blast against Sce and Rhto
blast = getBlast('hpo','hpo.faa',{'sce','rhto'},{'sce.faa','rhto_np11.faa'});
 
% Reconstruct model from both template at the same time, prefer sce reactions
modelHomology = getModelFromHomology({modelSce,modelRhto},blast,'hpo',{'sce','rhto'});

% Contract model, to merge identical reactions with different grRules
modelHomology = contractModel(modelHomology);

% Import template biomass equation to use as objective for later gapfilling
modelHomology = addRxnsGenesMets(modelHomology,modelSce,'r_4041');

% Only add exchange reactions for media components
mediumComps = {'r_1654','r_1672','r_1808','r_1832','r_1861','r_1992',...
   'r_2005','r_2060','r_2100','r_2111'};
modelHomology = addRxnsGenesMets(modelHomology,modelSce,mediumComps);

% Save workspace
save('HOMOLOGY.mat')
% load('HOMOLOGY.mat')

%% 3.4 PERFORM GAP-FILLING

% Use biomass production as obj func for gapfilling
modelHomology = setParam(modelHomology,'obj','r_4041',1);

% Set biomass production to arbitrary low flux, to force gap-filling to
% produce biomass.
modelHomology = setParam(modelHomology,'lb','r_4041',0.01);

% Set glycerol uptake at a higher value, to make sure that fluxes are high
% enough during gap-filling that they won't be ignored due to the tolerance
% of the MILP solver
modelHomology = setParam(modelHomology,'lb','r_1808',-10);

% From the Sce and Rhto models, remove all exchange reactions (the
% necessary ones we already added, don't want to add new ones)
modelSce2 = removeReactions(modelSce,getExchangeRxns(modelSce,'both'),true,true,true);
modelRhto2 = removeReactions(modelRhto,getExchangeRxns(modelRhto,'both'),true,true,true);

% Also, find which reactions in Rhto are already present in Sce (according
% to reaction ID, as Rhto model is based on Sce). No need to have duplicate
% reactions as suggestions during gap-filling, this reduces the solution
% space and helps to solve the MILP part of gap-filling.
rxns = intersect(modelSce2.rxns,modelRhto2.rxns);
modelRhto2 = removeReactions(modelRhto2,rxns);
 
% Run fillGaps function
[~,~, addedRxns, modelGapFilled, exitFlag]=fillGaps(modelHomology,{modelSce2,modelRhto2},false,true);

% Verify that model can now grow
sol=solveLP(modelGapFilled,1)
printFluxes(modelGapFilled,sol.x)

% Set maximum uptake of carbon source back to 1 mmol/gDCW*hr
modelHomology = setParam(modelHomology,'lb','r_1808',-1);

% Save workspace
save('GAPFILLING.mat')
% load('GAPFILLING.mat')

%% 3.5 CURATE BIOMASS STOICHIOMETRY

% Please refer to the text in order to obtain target-species-specific
% biomass stoichiometric coefficients for relevant fractions.

%% 3.5.1 DNA NUCLEOTIDES

rxnIndex= getIndexes(modelGapFilled,{'r_4050'},'rxns')
modelGapFilled.S(:,rxnIndex) % Use this line of code to determine the row number to edit for each nucleotide

modelGapFilled.S(338,rxnIndex)= -0.001893797; % Please note that the indexes are hardcoded and will not be the same in your case
modelGapFilled.S(343,rxnIndex)= -0.001738376;
modelGapFilled.S(354,rxnIndex)= -0.001738376;
modelGapFilled.S(373,rxnIndex)= -0.001893797;

%% 3.5.2 RNA NUCLEOTIDES

rxnIndex= getIndexes(modelGapFilled,{'r_4049'},'rxns')
modelGapFilled.S(:,rxnIndex)

modelGapFilled.S(230,rxnIndex)= -0.042547034;
modelGapFilled.S(295,rxnIndex)= -0.04503439;
modelGapFilled.S(433,rxnIndex)= -0.043433554;
modelGapFilled.S(904,rxnIndex)= -0.049036053;

%% 3.5.3 AMINO ACIDS

%{

0.089388158	Ala	A
0.067988923	Arg	R
0.062958735	Asn	N
0.079320115	Asp	D
0.017983071	Cys	C
0.055103267	Gln	Q
0.092543673	Glu	E
0.07422447	Gly	G
0.029418052	His	H
0.078633108	Ile	I
0.139071441	Leu	L
0.090151827	Lys	K
0.02898108	Met	M
0.060299163	Phe	F
0.062745261	Pro	P
0.112864926	Ser	S
0.074255724	Thr	T
0.015177252	Trp	W
0.04818661	Tyr	Y
0.087888537	Val	V

%}

rxnIndex= getIndexes(modelGapFilled,{'r_4047'},'rxns')
modelGapFilled.S(:,rxnIndex)

% Substrates (loaded tRNAs)
modelGapFilled.S(216,rxnIndex)= -0.089388158;
modelGapFilled.S(235,rxnIndex)= -0.067988923;
modelGapFilled.S(237,rxnIndex)= -0.062958735;
modelGapFilled.S(239,rxnIndex)= -0.079320115;
modelGapFilled.S(305,rxnIndex)= -0.017983071;
modelGapFilled.S(414,rxnIndex)= -0.055103267;
modelGapFilled.S(415,rxnIndex)= -0.092543673;
modelGapFilled.S(421,rxnIndex)= -0.07422447;
modelGapFilled.S(471,rxnIndex)= -0.029418052;
modelGapFilled.S(483,rxnIndex)= -0.078633108;
modelGapFilled.S(658,rxnIndex)= -0.139071441;
modelGapFilled.S(666,rxnIndex)= -0.090151827;
modelGapFilled.S(706,rxnIndex)= -0.02898108;
modelGapFilled.S(781,rxnIndex)= -0.060299163;
modelGapFilled.S(806,rxnIndex)= -0.062745261;
modelGapFilled.S(843,rxnIndex)= -0.112864926;
modelGapFilled.S(877,rxnIndex)= -0.074255724;
modelGapFilled.S(890,rxnIndex)= -0.015177252;
modelGapFilled.S(894,rxnIndex)= -0.04818661;
modelGapFilled.S(919,rxnIndex)= -0.087888537;

% Products (unloaded tRNAs should have same stoichiometry as corresponding
% substrate)
modelGapFilled.S(929,rxnIndex)= 0.089388158;
modelGapFilled.S(930,rxnIndex)= 0.067988923;
modelGapFilled.S(932,rxnIndex)= 0.062958735;
modelGapFilled.S(934,rxnIndex)= 0.079320115;
modelGapFilled.S(936,rxnIndex)= 0.017983071;
modelGapFilled.S(937,rxnIndex)= 0.055103267;
modelGapFilled.S(938,rxnIndex)= 0.092543673;
modelGapFilled.S(940,rxnIndex)= 0.07422447;
modelGapFilled.S(941,rxnIndex)= 0.029418052;
modelGapFilled.S(943,rxnIndex)= 0.078633108;
modelGapFilled.S(944,rxnIndex)= 0.139071441;
modelGapFilled.S(946,rxnIndex)= 0.090151827;
modelGapFilled.S(947,rxnIndex)= 0.02898108;
modelGapFilled.S(949,rxnIndex)= 0.060299163;
modelGapFilled.S(951,rxnIndex)= 0.062745261;
modelGapFilled.S(952,rxnIndex)= 0.112864926;
modelGapFilled.S(953,rxnIndex)= 0.074255724;
modelGapFilled.S(955,rxnIndex)= 0.015177252;
modelGapFilled.S(957,rxnIndex)= 0.04818661;
modelGapFilled.S(959,rxnIndex)= 0.087888537;


%% 3.5.4 LIPIDS

clear rxnIndex

exportToExcelFormat(modelGapFilled,'modelGFC.xlsx'); 

% Export model as a .xml file
exportModel(modelGapFilled, 'modelGapFilled.xml');

% Save workspace
save('BIOMASS.mat')
% load('BIOMASS.mat')

%% 3.6 SAVE TO GITHUB

% https://github.com/SysBioChalmers/Hansenula_polymorpha-GEM

% Create github repository similar to one linked above

% Upload first homologyModel draft

%% 3.7 SIMULATIONS


