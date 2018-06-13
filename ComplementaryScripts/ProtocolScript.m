%% PROTOCOL FOR THE RECONSTRUCTION OF A.FUMIGATUS GEM USING THE RAVEN TOOLBOX
%  AUTHOR: FRANCISCO ZORRILLA 

%{

Note: Some output from key functions will be displayed through comments,
be aware that this is done for didactic purposes. Its ok if slightly 
different outputs are observed when running the commands yourself,
especially if you use different protein FASTA files or if you use a
different GEM as a template for your desired model organism.

%}

%get new devel branch, fixed import model/miriams/chebi
%% Set working directory and load workspace (delete this section later)

cd C:\Users\zorrilla\Desktop\GEM

%% 1. RAVEN INSTALLATION
%     Dependencies: COBRA TOOLBOX (or libSBML), Gurobi, Gitwrapper

%{

Refer to https://github.com/SysBioChalmers/RAVEN/wiki/Installation for
details regarding installation process. Use the pathtool function to add
RAVEN, COBRA, and Gurobi subfolders to MATLAB path. ADD INFO ABOUT
GITWRAPPER, MEMOTE, TRAVIS(??)

%}

% Run the following lines to ensure that the installation was succesful:
global CBT_LP_SOLVER
CBT_LP_SOLVER = 'glpk';
checkInstallation;

%{

NOTE: Some MATLAB toolboxes may cause conflicitng errors with parsing excel
files. If checkInstallation returns "Checking if it is possible to parse a
model in Microsoft Excel format... FAILED", then please uninstall the
Text Analytics Toolbox.

%}

%% 2. PREPARATION

%{ 

Create a new folder to work in and store the following files.

Obtain high quality homologous organism GEM in .xml format. In this
protocol we will use the GEM for P. chrysogenum, which can be obtained from 
the supplementary materials of the following publication
http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002980.

Obtain protein FASTA of homologous organism. In this case we get the 
P. chrysogenum protein FASTA from the following link
ftp://ftp.ensemblgenomes.org/pub/fungi/release-39/fasta/fungi_ascomycota1_collection/penicillium_rubens_wisconsin_54_1255/pep/
using command line. Gunzip the downloaded file and rename the unzipped file
as pch_prot.

Obtain protein FASTA of target organism (A. fumigatus in this case)from NCBI at
ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/655/GCF_000002655.1_ASM265v1/GCF_000002655.1_ASM265v1_protein.faa.gz
using command line. Gunzip the downloaded file and rename the unzipped .faa file
as afm_prot.

NOTE: It is fine if the FASTA files have different extensions such as .fa,
.faa, .fasta, etc, as long as they are protein fasta. However, the gene
identifiers of the protein FASTA file of the homologous organism (in this
case P.chrysogenum) need to match the gene identifiers in the template
model. Check this by opening the fasta file using any text editor, and
comparing this with the template model (either in MATLAB after loading it
or simply by looking at the excel file of the model). In our case the FASTA
file we downloaded contains some extra information before/after the gene
identifiers, which will cause problems with the getModelFromHomology()
function. To avoid this, we need to remove these unnecessary annotations
using a text editor or command line. In our case we simply used Notepad++
and regular expressions to replace "(>.+)(Pc[0-9]{2}g[0-9]{5})(.*)" with ">\2"
to obtain a cleaned up FASTA file.

 %}
 
% Load P. chrysogenum GEM using importModel()
modelPch=importModel('iAl1006 v1.00.xml');

%{ 

modelPch = 

  struct with fields:

                id: 'iAL1006'
       description: 'Penicillium chrysogenum genome-scale model'
        annotation: [1×1 struct]
              rxns: {1632×1 cell}
              mets: {1235×1 cell}
                 S: [1235×1632 double]
                lb: [1632×1 double]
                ub: [1632×1 double]
               rev: [1632×1 double]
                 c: [1632×1 double]
                 b: [1235×1 double]
             comps: {5×1 cell}
         compNames: {5×1 cell}
       compOutside: {5×1 cell}
       compMiriams: {5×1 cell}
          rxnNames: {1632×1 cell}
           grRules: {1632×1 cell}
        rxnGeneMat: [1632×1006 double]
        subSystems: {1632×1 cell}
           eccodes: {1632×1 cell}
             genes: {1006×1 cell}
         geneComps: [1006×1 double]
       geneMiriams: {1006×1 cell}
    geneShortNames: {1006×1 cell}
          metNames: {1235×1 cell}
          metComps: [1235×1 double]
            inchis: {1235×1 cell}
       metFormulas: {1235×1 cell}
        metMiriams: {1235×1 cell}
%}

% Change model ID to "pch" for clarity
modelPch.id = 'pch';

% Set the objective function of the model to growth using setParam() and 
% confirm that the model is functional by using solveLP(). Reaction 1348
% correponds to growth in this particular model, this can be double 
% checked in the excel file.
modelPch=setParam(modelPch,'obj','r1348',1)

%{ 

modelPch = 

  struct with fields:

                id: 'pch'
       description: 'Penicillium chrysogenum genome-scale model'
        annotation: [1×1 struct]
              rxns: {1632×1 cell}
              mets: {1235×1 cell}
                 S: [1235×1632 double]
                lb: [1632×1 double]
                ub: [1632×1 double]
               rev: [1632×1 double]
                 c: [1632×1 double]
                 b: [1235×1 double]
             comps: {5×1 cell}
         compNames: {5×1 cell}
       compOutside: {5×1 cell}
       compMiriams: {5×1 cell}
          rxnNames: {1632×1 cell}
           grRules: {1632×1 cell}
        rxnGeneMat: [1632×1006 double]
        subSystems: {1632×1 cell}
           eccodes: {1632×1 cell}
             genes: {1006×1 cell}
         geneComps: [1006×1 double]
       geneMiriams: {1006×1 cell}
    geneShortNames: {1006×1 cell}
          metNames: {1235×1 cell}
          metComps: [1235×1 double]
            inchis: {1235×1 cell}
       metFormulas: {1235×1 cell}
        metMiriams: {1235×1 cell}

%}

% Next we need to set more restrictive boundaries for the exchange rate of 
% carbon sources to get a realistic growth rate. First, obtain exchange 
% reaction IDs using the getExchangeRxns()  function.
idx=getExchangeRxns(modelPch,'in');

% Constrain exchange reaction bounds using setParam() function
modelPch = setParam(modelPch,'eq',idx,0);

% Allow for uptake of some necessary metabolites
modelPch = setParam(modelPch,'ub',{'piIN','nh3IN','hno3IN','hno2IN',...
    'o2IN','slfIN','h2sIN','sulfurIN','h2so3IN','thmIN','pimIN'},1000);

% Allow for uptake of d-glucose as a carbon source
modelPch = setParam(modelPch,'ub',{'dglcIN'},1);

% Perform Flux Balance Analysis (FBA) using solveLP()
sol=solveLP(modelPch,1)

%{

sol = 

  struct with fields:

       x: [1632×1 double]
       f: -0.0857
    stat: 1
     msg: 'Optimal solution found'

%}

% Print exchange fluxes
printFluxes(modelPch,sol.x)

%{

EXCHANGE FLUXES:
c4odOUT	(production of 2-oxy-but-3-enoate):	8.5652e-06
bmOUT	(production of biomass):	0.085652
co2OUT	(production of CO2):	3.0316
h2oOUT	(production of H2O):	5.5818
dglcIN	(uptake of D-glucose):	1
piIN	(uptake of phosphate):	0.028169
nh3IN	(uptake of NH3):	0.59979
o2IN	(uptake of O2):	2.9526
h2sIN	(uptake of sulfide):	0.023108
sulfurIN	(uptake of sulfur):	8.5652e-06
thmIN	(uptake of thiamin):	8.5652e-06
pimIN	(uptake of pimelate):	8.5652e-06

%}

% Export P.chrysogenum model to excel for inspection
exportToExcelFormat(modelPch,'modelPch.xlsx');

%{ 

Open excel file in working folder to examine the model. The excel file
should have five different sheets named RXNS, METS, COMPS, GENES, and
MODEL.

%}

%% 3. CONSTRUCT DRAFT MODEL FROM HOMOLOGY, add second homology model from databases/fakeblast struct?

%{

To construct a draft model from homology we can use the protein FASTA 
files obtained in the PREPARATION section. Make sure gene identifiers match
as mentioned previously!

First, we perform a biderectional protein BLAST search using the getBlast()
function, to determine homologous proteins/enzyme catalyzed reactions in 
A.fumigatus and P.chrysogenum. This step results in a BLAST strucutre.


%}

%afmBlastStructure = getBlast('afm','afm_prot.faa',{'pch'},{'pch_prot.fa'});

%{

BLASTing "pch" against "afm"..
BLASTing "afm" against "pch"..

%}

% Save the BLAST structure to avoid unnecessarily re-running getBlast() in
% a future session
%save('afmBlastStructure.mat','afmBlastStructure');

% Next time use load() to access the BLAST structure
load('afmBlastStructure.mat');

% To avoid keeping unneccesary old genes, the models should not have
% "OR" relations in their grRules. To remove these use expandModel(),
% which separates reactions that were annotated with isoenzymes into
% multiple reactions. Note how there are more reactions in the expanded
% model than in the original.
modelPchExpanded=expandModel(modelPch)

%{

WARNING: Reaction r0123 contains nested and/or-relations. Large risk of errors

WARNING: Reaction r0409 contains nested and/or-relations. Large risk of errors

WARNING: Reaction r0465 contains nested and/or-relations. Large risk of errors

WARNING: Reaction r0631 contains nested and/or-relations. Large risk of errors

WARNING: Reaction r0651 contains nested and/or-relations. Large risk of errors

WARNING: Reaction r0673 contains nested and/or-relations. Large risk of errors

WARNING: Reaction r0702 contains nested and/or-relations. Large risk of errors

WARNING: Reaction r0705 contains nested and/or-relations. Large risk of errors


modelPchExpanded = 

  struct with fields:

                id: "pch"
       description: 'Penicillium chrysogenum genome-scale model'
        annotation: [1×1 struct]
              rxns: {2715×1 cell}
              mets: {1235×1 cell}
                 S: [1235×2715 double]
                lb: [2715×1 double]
                ub: [2715×1 double]
               rev: [2715×1 double]
                 c: [2715×1 double]
                 b: [1235×1 double]
             comps: {5×1 cell}
         compNames: {5×1 cell}
       compOutside: {5×1 cell}
       compMiriams: {5×1 cell}
          rxnNames: {2715×1 cell}
           grRules: {2715×1 cell}
        rxnGeneMat: [2715×1006 double]
        subSystems: {2715×1 cell}
           eccodes: {2715×1 cell}
             genes: {1006×1 cell}
         geneComps: [1006×1 double]
       geneMiriams: {1006×1 cell}
    geneShortNames: {1006×1 cell}
          metNames: {1235×1 cell}
          metComps: [1235×1 double]
            inchis: {1235×1 cell}
       metFormulas: {1235×1 cell}
        metMiriams: {1235×1 cell}

%}

% Obtain model for A.fumigatus using BLAST structure and the expanded
% P.chrysogenum model using getModelFromHomology().
modelAfm=getModelFromHomology({modelPchExpanded},afmBlastStructure,'afm',{},1,false,10^-20,100,35)

%{

modelAfm = 

  struct with fields:

                     id: 'afm'
            description: 'Generated by getModelFromHomology using pch'
             annotation: [1×1 struct]
                   rxns: {2191×1 cell}
                   mets: {1078×1 cell}
                      S: [1078×2191 double]
                     lb: [2191×1 double]
                     ub: [2191×1 double]
                    rev: [2191×1 double]
                      c: [2191×1 double]
                      b: [1078×1 double]
                  comps: {5×1 cell}
              compNames: {5×1 cell}
            compOutside: {5×1 cell}
            compMiriams: {5×1 cell}
               rxnNames: {2191×1 cell}
                grRules: {2191×1 cell}
             rxnGeneMat: [2191×1058 double]
             subSystems: {2191×1 cell}
                eccodes: {2191×1 cell}
                  genes: {1058×1 cell}
              geneComps: [1058×1 double]
            geneMiriams: {966×1 cell}
               metNames: {1078×1 cell}
               metComps: [1078×1 double]
                 inchis: {1078×1 cell}
            metFormulas: {1078×1 cell}
             metMiriams: {1078×1 cell}
               rxnNotes: {2191×1 cell}
    rxnConfidenceScores: [2191×1 double]

%}

% Contract model to combine the reactions were seperated by expandModel() 
% above. Isoenzymes are now annotated to the same reaction, with 'OR' gene 
% relationships. Note how this model has less reactions than the one
% obtained using the getModelFromHomology() function.
modelAfm=contractModel(modelAfm)

%{

modelAfm = 

  struct with fields:

                     id: 'afm'
            description: 'Generated by getModelFromHomology using pch'
             annotation: [1×1 struct]
                   rxns: {1149×1 cell}
                   mets: {1078×1 cell}
                      S: [1078×1149 double]
                     lb: [1149×1 double]
                     ub: [1149×1 double]
                    rev: [1149×1 double]
                      c: [1149×1 double]
                      b: [1078×1 double]
                  comps: {5×1 cell}
              compNames: {5×1 cell}
            compOutside: {5×1 cell}
            compMiriams: {5×1 cell}
               rxnNames: {1149×1 cell}
                grRules: {1149×1 cell}
             rxnGeneMat: [1149×1058 double]
             subSystems: {1149×1 cell}
                eccodes: {1149×1 cell}
                  genes: {1058×1 cell}
              geneComps: [1058×1 double]
            geneMiriams: {966×1 cell}
               metNames: {1078×1 cell}
               metComps: [1078×1 double]
                 inchis: {1078×1 cell}
            metFormulas: {1078×1 cell}
             metMiriams: {1078×1 cell}
               rxnNotes: {1149×1 cell}
    rxnConfidenceScores: [1149×1 double]

%}

% Remove reactions with EXP_ in the reaction IDs. These were introduced 
% by expandModel() to allow duplicate reactions when isoenzymes were 
% annotated. However, we contracted the model here again, so any remaining
% reactions with EXP_ in the reaction ID is an artifact (where an ortholog 
% was found for one of the isoenzymes, but not for the others). Note that
% in this case there were no artifact reactions to remove.
modelAfm.rxns=regexprep(modelAfm.rxns,'_EXP_.','') 

%{

modelAfm = 

  struct with fields:

                     id: 'afm'
            description: 'Generated by getModelFromHomology using pch'
             annotation: [1×1 struct]
                   rxns: {1149×1 cell}
                   mets: {1078×1 cell}
                      S: [1078×1149 double]
                     lb: [1149×1 double]
                     ub: [1149×1 double]
                    rev: [1149×1 double]
                      c: [1149×1 double]
                      b: [1078×1 double]
                  comps: {5×1 cell}
              compNames: {5×1 cell}
            compOutside: {5×1 cell}
            compMiriams: {5×1 cell}
               rxnNames: {1149×1 cell}
                grRules: {1149×1 cell}
             rxnGeneMat: [1149×1058 double]
             subSystems: {1149×1 cell}
                eccodes: {1149×1 cell}
                  genes: {1058×1 cell}
              geneComps: [1058×1 double]
            geneMiriams: {966×1 cell}
               metNames: {1078×1 cell}
               metComps: [1078×1 double]
                 inchis: {1078×1 cell}
            metFormulas: {1078×1 cell}
             metMiriams: {1078×1 cell}
               rxnNotes: {1149×1 cell}
    rxnConfidenceScores: [1149×1 double]

%}

% Add all exchange reactions. These were not gene annotated, and therefore 
% not added by the getModelFromHomology() function. We might not require 
% all exchange rxns, but we can easily remove unconnected ones later.
exRxns=getExchangeRxns(modelPch);                                                %% NOTE MAYBE BETTER TO DO THIS AND
modelAfm=addRxnsGenesMets(modelAfm,modelPch,exRxns,false,'Modeling reaction',1); %% REST OF THIS SECTION AFTER MERGING

%{

The following metabolites were already present in the model and will not be added:
aminoacetaldehyde[c]
glycine betaine[c]
nicotinate[c]
1,3-beta-D-glucan[e]
acetate[e]
adenine[e]
alpha,alpha-trehalose[e]
alpha-D-glucose[e]
beta-D-glucose[e]
cellobiono-1,5-lactone[e]
cellobiose[e]
cellulose[e]
chitin[e]
chitobiose[e]
chitosan[e]
choline[e]
cytosine[e]
D-fructose[e]
D-galactose[e]
D-glucosamine[e]
D-glucose[e]
D-mannose[e]
D-xylose[e]
gamma-aminobutyrate[e]
glycine[e]
glycogen[e]
guanine[e]
H2O[e]
H2O2[e]
lactose[e]
L-alanine[e]
L-arginine[e]
L-asparagine[e]
L-aspartate[e]
L-cysteine[e]
L-dopaquinone[e]
L-glutamate[e]
L-glutamine[e]
L-histidine[e]
L-isoleucine[e]
L-leucine[e]
L-lysine[e]
L-methionine[e]
L-phenylalanine[e]
L-proline[e]
L-serine[e]
L-threonine[e]
L-tryptophan[e]
L-tyrosine[e]
L-valine[e]
maltose[e]
mannans[e]
melibiose[e]
myo-inositol[e]
N-acetyl-D-glucosamine[e]
NH3[e]
nitrate[e]
O2[e]
raffinose[e]
starch[e]
sucrose[e]
sulfate[e]
thiamin[e]
uracil[e]
uridine[e]
xylans[e]

The following metabolites will be added to the model:
N-[L-5-amino-5-carboxypentanoyl]-L-cysteinyl-D-valine[e]
(S)-lactate[e]
(S)-malate[e]
2-hydroxybenzylpenicillin[e]
2-oxoglutarate[e]
2-oxy-but-3-enoate[e]
4-aminobenzoate[e]
4-hydroxyphenoxyacetate[e]
4-hydroxyphenoxymethylpenicillin[e]
6-aminopenicillanate[e]
6-oxopiperidine-2-carboxylate[e]
8-hydroxypenillic acid[e]
alpha-L-arabinan[e]
anthranilate[e]
artificial penicillin[e]
artificial protein[e]
benzylpenicillin[e]
benzylpenicilloic acid[e]
beta-alanine[e]
biomass[e]
butyrate[e]
citrate[e]
CO2[e]
cyanate[e]
D-arabinitol[e]
D-arabinose[e]
decanoate[e]
D-gluconate[e]
D-glucono-1,5-lactone[e]
D-mannitol[e]
D-ribose[e]
ethanol[e]
ethylnitronate[e]
FMN[e]
formate[e]
fumarate[e]
glutathione[e]
glycerol[e]
glycolate[e]
H2S[e]
heptadecanoate[e]
heptadecenoate[e]
heptanoate[e]
hexadecadienoate[e]
hexadecenoate[e]
hexanoate[e]
hexanoylpenicillin[e]
hypoxanthine[e]
icosanoate[e]
isocitrate[e]
isopenicillin N[e]
L-2-aminoadipate[e]
L-arabinitol[e]
L-arabinose[e]
laurate[e]
L-citrulline[e]
L-cystine[e]
L-homocysteine[e]
L-iditol[e]
L-ornithine[e]
L-ribulose[e]
L-sorbose[e]
maltotriose[e]
methanol[e]
myristate[e]
nicotinamide[e]
nitrite[e]
nonanoate[e]
octadecadienoate[e]
octadecatrienoate[e]
octadecenoate[e]
octanoate[e]
octanoylpenicillin[e]
oxalate[e]
oxaloacetate[e]
palmitate[e]
pentadecanoate[e]
phenoxyacetate[e]
phenoxymethylpenicillin[e]
phenoxymethylpenicilloic acid[e]
phenylacetate[e]
phosphate[e]
pimelate[e]
propionate[e]
pyruvate[e]
quinolinate[e]
stearate[e]
succinate[e]
sulfite[e]
sulfur[e]
urea[e]
valerate[e]
xanthine[e]
xylitol[e]

Number of metabolites added to the model:
94


Number of reactions added to the model:
161

%}

% Same as exchange reactions, add all non-gene annotated transport reactions
noGeneIdx=find(cellfun(@isempty,modelPch.grRules)); % Find rxns with no genes
rxnIdx=regexp(modelPch.rxnNames,'(transport)|(diffusion)|(carrier)|(shuttle)'); % Find rxnNames containing 'transport' or 'diffusion'
rxnIdx=find(~cellfun('isempty',rxnIdx)); % Find their indices
rxnIdx=intersect(rxnIdx,noGeneIdx); % Keep the ones without gene notation
exRxns=modelPch.rxns(rxnIdx); % Obtain reaction IDs
modelAfm=addRxnsGenesMets(modelAfm,modelPch,exRxns,false,...
    'Modeling reaction required for intercellular transport, gene unknown',1); % Finally add reactions and metabolites

%{

The following metabolites were already present in the model and will not be added:
(S)-lactate[c]
(S)-malate[c]
2-hydroxyphenylacetate[c]
4-hydroxyphenoxyacetate[c]
acetaldehyde[c]
acetate[c]
acetoacetate[c]
acetylcarnitine[c]
acyl-[acp][c]
AMP[c]
ATP[c]
butyrate[c]
carnitine[c]
CO2[c]
decanoate[c]
diphosphate[c]
ethanol[c]
formate[c]
H2O[c]
H2S[c]
heptadecanoate[c]
heptadecenoate[c]
heptanoate[c]
hexadecadienoate[c]
hexadecenoate[c]
hexadecenoyl-[acp][c]
hexanoate[c]
icosanoate[c]
icosanoyl-[acp][c]
isocitrate[c]
laurate[c]
methanol[c]
myristate[c]
N-[L-5-amino-5-carboxypentanoyl]-L-cysteinyl-D-valine[e]
NH3[c]
nonanoate[c]
O2[c]
octadecadienoate[c]
octadecadienoyl-[acp][c]
octadecatrienoate[c]
octadecenoate[c]
octadecenoyl-[acp][c]
octanoate[c]
palmitate[c]
palmitoyl-[acp][c]
pentadecanoate[c]
pentadecanoyl-[acp][c]
phenoxyacetate[c]
phenylacetate[c]
stearate[c]
stearoyl-[acp][c]
(S)-lactate[e]
2-oxy-but-3-enoate[e]
4-hydroxyphenoxyacetate[e]
acetate[e]
butyrate[e]
CO2[e]
decanoate[e]
ethanol[e]
formate[e]
H2O[e]
heptadecanoate[e]
heptadecenoate[e]
heptanoate[e]
hexadecadienoate[e]
hexadecenoate[e]
hexanoate[e]
icosanoate[e]
laurate[e]
methanol[e]
myristate[e]
nonanoate[e]
O2[e]
octadecadienoate[e]
octadecatrienoate[e]
octadecenoate[e]
octanoate[e]
palmitate[e]
pentadecanoate[e]
phenoxyacetate[e]
phenylacetate[e]
stearate[e]
(S)-malate[m]
acetaldehyde[m]
acetate[m]
acetoacetate[m]
CO2[m]
ethanol[m]
H2O[m]
H2S[m]
isocitrate[m]
NH3[m]
O2[m]
2-hydroxyphenylacetate[p]
4-hydroxyphenoxyacetate[p]
AMP[p]
ATP[p]
decanoate[p]
diphosphate[p]
H2O[p]
heptadecanoate[p]
icosanoate[p]
laurate[p]
myristate[p]
O2[p]
octanoate[p]
palmitate[p]
pentadecanoate[p]
phenoxyacetate[p]
phenylacetate[p]
stearate[p]

The following metabolites will be added to the model:
2-oxy-but-3-enoate[c]
N-[L-5-amino-5-carboxypentanoyl]-L-cysteinyl-D-valine[c]
acetylcarnitine[p]
carnitine[p]
CO2[p]

Number of metabolites added to the model:
5

Number of reactions added to the model:
61
 
%}

% Ensure biomass formation and transport reactions are present in model
modelAfm = addRxnsGenesMets(modelAfm, modelPch,{'r1463','r1348'})

% Remove unnecessary geneMiriams field in modelAfm strucutre. This will be
% done automatically in future update of getModelFromHomology() function.
% Rename model to highlight that it was obtained from homology.
modelAfmHomology = rmfield(modelAfm,'geneMiriams');
    
% Export model to excel for inspection
exportToExcelFormat(modelAfmHomology,'modelAfmHomology.xlsx'); 

% Add spontaneous bicarbonate formation DOUBLE CHECK IF THIS APPLIES TO
% P.CH MODEL. Looks like bicarbonate formation is not present in Pch model,
% maybe good to include in final draft of A.fumigatus model???????????
%exRxns={'r_1664','r_1665','r_1667','r_1668'};
%modelCal2=addRxnsGenesMets(modelCal2,modelSce,exRxns,true,...
%    'Manual curation','1'); % Add reactions and metabolites

% Clean up workspace
clear rxnIdx noGeneIdx exRxns idx modelPchExpanded afmBlastStructure sol modelAfm

% Check that model can grow
modelAfmHomology = setParam(modelAfmHomology,'obj','r1348',1);
sol=solveLP(modelAfmHomology,1);
printFluxes(modelAfmHomology,sol.x);
% probably expected that model will be unable to grow

%% 4. GITHUB (most of this section will be done outside of MATLAB, meet with benjamin to discuss)

% https://github.com/SysBioChalmers/Aspergillus_fumigatus-GEM

%{

%questions: 
%   what should be included in repo? look at S.ce repo
%   what is the git wrapper? add to dependencies
%   how should github be practiaclly used throughout the process of reconstruction (after
%   running each section? after reconstruction using different methods?)
%   use memote-travis

% make repo manually 
% git wrapper
% meet w/benjamin on friday @ 13:30
% look at S.coelicolor repo, for commit script

%}

% Install wrapper

% Create github repository similar to one linked above

% Upload first homologyModel draft

%% 5. CONSTRUCT DRAFT MODELS FROM KEGG - check excel sheets, look at subsystem

%{

In this section we make use of the KEGG-based resconstruction approaches
offered by the RAVEN toolbox. The first approach is based on hidden markov 
models (HHMs), while the second approach is based on annotation.

%}

% Create model using KEGG & HMM
modelAfmKeggHMM = getKEGGModelForOrganism('afm','afm_prot.faa','euk100_kegg82','output',false,false,false,10^-30,0.8,0.3,-1)

%{

Downloading HMMs archive file...
Extracting HMMs archive file...
NOTE: Importing KEGG reactions from C:/Users/zorrilla/Desktop/GEM/RAVEN-master/external/kegg/keggRxns.mat.
KEGG reactions loaded
NOTE: Importing KEGG genes from C:/Users/zorrilla/Desktop/GEM/RAVEN-master/external/kegg/keggGenes.mat.
KEGG genes loaded
NOTE: Importing KEGG metabolites from C:/Users/zorrilla/Desktop/GEM/RAVEN-master/external/kegg/keggMets.mat.
KEGG metabolites loaded
Completed generation of KEGG model
NOTE: Importing KEGG phylogenetic distance matrix from C:/Users/zorrilla/Desktop/GEM/RAVEN-master/external/kegg/keggPhylDist.mat.
Completed creation of phylogenetic distance matrix
Completed generation of multi-FASTA files
Completed multiple alignment of sequences
Completed generation of HMMs
Completed matching to HMMs

modelAfmKeggHMM = 

  struct with fields:

             id: 'afm'
    description: 'Automatically generated from KEGG database'
           rxns: {1842×1 cell}
       rxnNames: {1842×1 cell}
        eccodes: {1842×1 cell}
     subSystems: {1842×1 cell}
     rxnMiriams: {1842×1 cell}
              S: [1880×1842 double]
           mets: {1880×1 cell}
            rev: [1842×1 double]
             ub: [1842×1 double]
             lb: [1842×1 double]
              c: [1842×1 double]
              b: [1880×1 double]
          genes: {1664×1 cell}
     rxnGeneMat: [1842×1664 double]
       metNames: {1880×1 cell}
    metFormulas: {1880×1 cell}
         inchis: {1880×1 cell}
     metMiriams: {1880×1 cell}
          comps: {'s'}
      compNames: {'System'}
    compOutside: {''}
       metComps: [1880×1 double]
        grRules: {1842×1 cell}

%}
        
exportToExcelFormat(modeAfmKeggHMM,'modelAfmKeggHMM.xlsx') % Inspect model

% Create model using KEGG & annotation
modelAfmKeggAnnot =getKEGGModelForOrganism('afm','','','',false,false,false)   %check if this worked corrctly

%{

NOTE: Importing KEGG reactions from C:/Users/zorrilla/Desktop/GEM/RAVEN-master/external/kegg/keggRxns.mat.
KEGG reactions loaded
NOTE: Importing KEGG genes from C:/Users/zorrilla/Desktop/GEM/RAVEN-master/external/kegg/keggGenes.mat.
KEGG genes loaded
NOTE: Importing KEGG metabolites from C:/Users/zorrilla/Desktop/GEM/RAVEN-master/external/kegg/keggMets.mat.
KEGG metabolites loaded
Completed generation of KEGG model

modelAfmKeggAnnot = 

  struct with fields:

             id: 'afm'
    description: 'Automatically generated from KEGG database'
           rxns: {1404×1 cell}
       rxnNames: {1404×1 cell}
        eccodes: {1404×1 cell}
     subSystems: {1404×1 cell}
     rxnMiriams: {1404×1 cell}
              S: [1505×1404 double]
           mets: {1505×1 cell}
            rev: [1404×1 double]
             ub: [1404×1 double]
             lb: [1404×1 double]
              c: [1404×1 double]
              b: [1505×1 double]
          genes: {1110×1 cell}
     rxnGeneMat: [1404×1110 double]
       metNames: {1505×1 cell}
    metFormulas: {1505×1 cell}
         inchis: {1505×1 cell}
     metMiriams: {1505×1 cell}
          comps: {'s'}
      compNames: {'System'}
    compOutside: {''}
       metComps: [1505×1 double]
        grRules: {1404×1 cell}
    geneMiriams: {1110×1 cell}



%}

exportToExcelFormat(modeAfmKeggAnnot,'modelAfmKeggAnnot.xlsx') % Inspect model

%% 6. CONSTRUCT DRAFT MODEL FROM METACYC --- REVIEW
%
modelAfmMC=getMetaCycModelForOrganism('afm','afm_prot.faa',false,false,false,100,55) %double check parameter w ed

%{

NOTE: Importing MetaCyc reactions from C:/Users/zorrilla/Desktop/GEM/RAVEN-master/external/metacyc/metaCycRxns.mat.
MetaCyc reactions loaded
NOTE: Importing MetaCyc enzymes and reaction-enzyme association from C:/Users/zorrilla/Desktop/GEM/RAVEN-master/external/metacyc/metaCycEnzymes.mat.
MetaCyc enzymes loaded
NOTE: Importing MetaCyc metabolites from C:/Users/zorrilla/Desktop/GEM/RAVEN-master/external/metacyc/metaCycMets.mat.
MetaCyc compounds loaded
The full MetaCyc model loaded
BLASTing "MetaCyc" against "afm"..
BLASTing "afm" against "MetaCyc"..
Completed searching against MetaCyc protein sequences.
NOTE: Importing MetaCyc metabolites from C:/Users/zorrilla/Desktop/GEM/RAVEN-master/external/metacyc/metaCycMets.mat.

modelAfmMC = 

  struct with fields:

                     id: 'afm'
            description: 'Generated by homology with MetaCyc database'
                   rxns: {1557×1 cell}
               rxnNames: {1557×1 cell}
                eccodes: {1557×1 cell}
             subSystems: {1557×1 cell}
             rxnMiriams: {1557×1 cell}
          rxnReferences: {1557×1 cell}
                     lb: [1557×1 double]
                     ub: [1557×1 double]
                    rev: [1557×1 double]
                      c: [1557×1 double]
              equations: {1557×1 cell}
                  genes: {1987×1 cell}
             rxnGeneMat: [1557×1987 double]
    rxnConfidenceScores: [1557×1 double]
                grRules: {1557×1 cell}
                      S: [2015×1557 double]
                   mets: {2015×1 cell}
               metNames: {2015×1 cell}
            metFormulas: {2015×1 cell}
             metCharges: [2015×1 double]
                      b: [2015×1 double]
                 inchis: {2015×1 cell}
             metMiriams: {2015×1 cell}
                  comps: {'s'}
              compNames: {'System'}
               metComps: [2015×1 double]

%}

exportToExcelFormat(modelAfmMC,'modelAfmMC.xlsx');

%% 7. COMBINE MODELS -rerun with new MC model
% 

%{

This is perhaps the trickiest part of the reconstruction procedure, as
merging models can be extremely tedious and time consuming if done in a 
manual fashion. On the other hand, taking an automated approach may result
in many errors when matching metabolites and/or reactions. One approach 
which has been developed to tackle this challange is the modelBorgifier 
Toolbox, which is part of the COBRA toolbox.

%}

% SHOULD I DO GAP FILLING BEFROE LOCALIZATION? WHAT MODEL DO I USE AS
% TEMPLATE? DO I NEED TO ADD EXCH RXNS TO MCKEGG MODEL BEFORE GAPFILL?

% Use mergeModels() function to combine the two KEGG-based approach models
modelAfmKegg= mergeModels({modelAfmKeggHMM,modelAfmKeggAnnot}) %these models have missing names for reactions

%{

modelAfmKegg = 

  struct with fields:

             id: 'MERGED'
    description: ''
           rxns: {3265×1 cell}
       rxnNames: {3265×1 cell}
        eccodes: {3265×1 cell}
     subSystems: {3265×1 cell}
     rxnMiriams: {3265×1 cell}
              S: [1924×3265 double]
           mets: {1924×1 cell}
            rev: [3265×1 double]
             ub: [3265×1 double]
             lb: [3265×1 double]
              c: [3265×1 double]
              b: [1924×1 double]
          genes: {2795×1 cell}
     rxnGeneMat: [3265×2795 double]
       metNames: {1924×1 cell}
    metFormulas: {1924×1 cell}
         inchis: {1924×1 cell}
     metMiriams: {1924×1 cell}
          comps: {'s'}
      compNames: {'System'}
    compOutside: {''}
       metComps: [1924×1 double]
        grRules: {3265×1 cell}
        rxnFrom: {3265×1 cell}
        metFrom: {1924×1 cell}
       geneFrom: {2795×1 cell}
    geneMiriams: {2795×1 cell}

%}

exportToExcelFormat(modelAfmKegg,'modelAfmKegg.xlsx');

%use combineMetaCycKEGGModels() instead, NOTE these dont have compartmens,
%use predictLocalization() (need linux) function for this (talk to simonas). Look at
%metanetX for matching/annotation (Hao), look at tutorial modelBorgefier script
%part of cobra toolbox(reconstruction/comparison)

modelAfmMCKegg2 = combineMetaCycKEGGModels(modelAfmMC55,modelAfmKegg)

%{

NOTE: A total of 1132 reactions in the KEGG model were mapped to MetaCyc.
NOTE: 626 reactions already in MetaCyc model, 506 will be combined.
NOTE: A shrinked KEGG model with 2133 reactions and 1764 metabolites was obtained.
NOTE: A total of 943 metabolites from the shrinked KEGG model were mapped to MetaCyc again.
NOTE: Importing KEGG metabolites from C:/Users/zorrilla/Desktop/GEM/RAVEN-master/external/kegg/keggMets.mat.

modelAfmMCKegg = 

  struct with fields:

                     id: 'COMBINED'
            description: 'Combined model from MetaCyc and KEGG draft models'
                   rxns: {4196×1 cell}
               rxnNames: {4196×1 cell}
                eccodes: {4196×1 cell}
             subSystems: {4196×1 cell}
             rxnMiriams: {4196×1 cell}
          rxnReferences: {4196×1 cell}
                     lb: [4196×1 double]
                     ub: [4196×1 double]
                    rev: [4196×1 double]
                      c: [4196×1 double]
              equations: {4196×1 cell}
                  genes: {3513×1 cell}
                grRules: {4196×1 cell}
                   mets: {3353×1 cell}
             metCharges: [3353×1 double]
                      b: [3353×1 double]
                 inchis: {3353×1 cell}
             metMiriams: {3353×1 cell}
                  comps: {'s'}
              compNames: {'System'}
               metComps: [3353×1 double]
                rxnFrom: {4196×1 cell}
               geneFrom: {3513×1 cell}
            grRulesKEGG: {4196×1 cell}
                      S: [3353×4196 double]
                metFrom: {3353×1 cell}
               metNames: {3353×1 cell}
            metFormulas: {3353×1 cell}
    rxnConfidenceScores: [4196×1 double]
             rxnGeneMat: [4196×3513 double]

%}

exportToExcelFormat(modelAfmMCKegg,'modelAfmMCKegg.xlsx');

% combine MCKegg model with Homology model (using modelBorgifier)

save(modelAfmMCKegg,'modelAfmMCKegg.mat')
save('modelAfmHomology.mat')

%% localization
 % remove localization in homology mergeCompartments() (keep unconstrained), merge models,gap fill ,predict localization
    %use homology/pch modelsz for gap filling 
    
    %do not want all of the compartments from cello, use copyToComps()
    %function to move rxns arounnd/then delete unwanted comps. keep cytosol,
    %mitochondria, peroxisome, extracellular; everything else is cytosolic
    
    %DO:
    %merge compartments in AfmHomology/Pch, fillGaps, predict localization
    %(compare to original AfmHomology to see if predict localization is working), 
    % then rerun this starting with AfmHomologyKeggMC
    
    %predict localization on KEGG, MC models, if works merge with Homology

testGSS = parseScores('seqcat.txt','cello');
[modelAfmHomologyMergeComps, deletedRxnsMergeComps, duplicateRxnsMergeComp]=mergeCompartments(modelAfmHomology,true);
modelAfmMerged = mergeModels({modelAfmHomologyMergeComps,modelAfmMCKegg});
[newConnected, cannotConnect, addedRxns, newModel, exitFlag]=fillGaps(modelAfmMerged,modelPch,false,true);
%^ try specifying objective function for biomass equation,
%modelAfmMerged/modelPch should both either be compartment merged or not.
[outModel, geneLocalization, transportStruct, score, removedRxns] = predictLocalization(modelAfmMerged,testGSS,'cytop') %doesnt work
%{
Matrix index is out of range for deletion.

Error in removeMets (line 101)
        reducedModel.metCharges(indexesToDelete)=[];

Error in removeReactions (line 114)
            reducedModel=removeMets(reducedModel,unUsedMets,false,false,false,removeUnusedComps);

Error in predictLocalization (line 167)
        model=removeReactions(model,I,true,true);
%}

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

