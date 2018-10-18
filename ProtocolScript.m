%% PROTOCOL FOR THE RECONSTRUCTION OF A.FUMIGATUS GEM USING THE RAVEN TOOLBOX
%  AUTHOR: FRANCISCO ZORRILLA 

%{

Note: Some output from key functions will be displayed through comments,
be aware that this is done for didactic purposes. Its ok if slightly 
different outputs are observed when running the commands yourself,
especially if you use different protein FASTA files or if you use a
different GEM as a template for your desired model organism.

%}


%% Set working directory and load workspace (delete this section later)

cd C:\Users\zorrilla\Desktop\RAVENprotocol

%% 1. RAVEN INSTALLATION
%     Dependencies: COBRA TOOLBOX (or libSBML), Gurobi, Gitwrapper?

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

% Perform Flux Balance Analysis using solveLP()
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

%{ 

i.   Create a new folder to work in and store the following files.

ii.  Obtain high quality homologous organism GEM in .xml format. In this
     protocol we will use the GEM for P. chrysogenum, which can be obtained from 
     the supplementary materials of the following publication
     http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002980.

iii. Obtain protein FASTA of homologous organism. In this case we get the 
     P. chrysogenum protein FASTA from the following link
     ftp://ftp.ensemblgenomes.org/pub/fungi/release-39/fasta/fungi_ascomycota1_collection/penicillium_rubens_wisconsin_54_1255/pep/Penicillium_rubens_wisconsin_54_1255.PenChr_Nov2007.pep.all.fa.gz     using command line. Gunzip the downloaded file and rename the unzipped file
     as pch_prot.

iv.  Obtain protein FASTA of target organism (A. fumigatus in this case)from NCBI at
     ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/655/GCF_000002655.1_ASM265v1/GCF_000002655.1_ASM265v1_protein.faa.gz
     using command line. Gunzip the downloaded file and rename the unzipped .faa file
     as afm_prot.

NOTE: It is fine if FASTA files have different extensions such as .fa,
.faa, .fasta, etc, as long as they are protein fasta files. However, the gene
identifiers of the protein FASTA file of the homologous organism (in this
case P.chrysogenum) need to match the gene identifiers in the template
model. Check this by opening the fasta file using any text editor, and
comparing this with the template model (either in MATLAB after loading it
or simply by looking at the excel file of the model). In our case the FASTA
file we downloaded contains some extra information before/after the gene
identifiers, which will cause problems with the getModelFromHomology()
function. To avoid this, we need to remove these unnecessary annotations
using a text editor or command line. In our case we simply used Notepad++
and regular expressions to find "(>.+)(Pc[0-9]{2}g[0-9]{5})(.*)" and replace 
with ">\2" to obtain a cleaned up FASTA file.

 %}

%% 3. CONSTRUCT DRAFT MODEL FROM HOMOLOGY, ASK ED ABT second homology model from databases/fakeblast struct? (metaphors)
          %CONSIDER MERGING COMPS OF TEMPLATE MODEL BEFORE GETFROMHOMOLOGY
          %-attempting- ask ED if worth trying
%{

To construct a draft model from homology we can use the protein FASTA 
files obtained in the PREPARATION section. Make sure gene identifiers match
as mentioned previously!

First, we perform a biderectional protein BLAST search using the getBlast()
function, to determine homologous proteins/enzyme catalyzed reactions in 
A.fumigatus and P.chrysogenum. This step results in a BLAST strucutre.


%}

afmBlastStructure = getBlast('afm','afm_prot.fa',{'pch'},{'pch_prot.fa'});

%{

BLASTing "pch" against "afm"..
BLASTing "afm" against "pch"..

%}

% Save the BLAST structure to avoid unnecessarily re-running getBlast() in
% a future session
save('afmBlastStructure.mat','afmBlastStructure');

% Next time use load() to access the BLAST structure
load('afmBlastStructure.mat');

% Merge compartments of template model to allow for easier model merging
% later on, localization will be re-established after merging
[modelPch, deletedRxns, duplicateRxns]=mergeCompartments(modelPch,true)


% To avoid keeping unneccesary old genes, the models should not have
% 'OR' relations in their grRules. To remove these use expandModel(),
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
%exRxns=getExchangeRxns(modelPch);                                                %% NOTE MAYBE BETTER TO DO THIS AND
%modelAfm=addRxnsGenesMets(modelAfm,modelPch,exRxns,false,'Modeling reaction',1); %% REST OF THIS SECTION AFTER MERGING add exch after gap filling

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
%noGeneIdx=find(cellfun(@isempty,modelPch.grRules)); % Find rxns with no genes
%rxnIdx=regexp(modelPch.rxnNames,'(transport)|(diffusion)|(carrier)|(shuttle)'); % Find rxnNames containing 'transport' or 'diffusion'
%rxnIdx=find(~cellfun('isempty',rxnIdx)); % Find their indices
%rxnIdx=intersect(rxnIdx,noGeneIdx); % Keep the ones without gene notation
%exRxns=modelPch.rxns(rxnIdx); % Obtain reaction IDs
%modelAfm=addRxnsGenesMets(modelAfm,modelPch,exRxns,false,...
%    'Modeling reaction required for intercellular transport, gene unknown',1); % Finally add reactions and metabolites

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

% Remove unnecessary geneMiriams field in modelAfm strucutre. This will be
% done automatically in future update of getModelFromHomology() function.
% Rename model to highlight that it was obtained from homology.
modelAfmHomology = rmfield(modelAfm,'geneMiriams');
%modelAfmHomology = rmfield(modelAfm,'geneComps');    

% Export model to excel for inspection
exportToExcelFormat(modelAfmHomology,'modelAfmHomology.xlsx'); 

% Add spontaneous bicarbonate formation DOUBLE CHECK IF THIS APPLIES TO
% P.CH MODEL. Looks like bicarbonate formation is not present in Pch model,
% maybe good to include in final draft of A.fumigatus model
%exRxns={'r_1664','r_1665','r_1667','r_1668'};
%modelCal2=addRxnsGenesMets(modelCal2,modelSce,exRxns,true,...
%    'Manual curation','1'); % Add reactions and metabolites

% Clean up workspace
clear rxnIdx noGeneIdx exRxns idx modelPchExpanded afmBlastStructure sol modelAfm


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

%% 5. CONSTRUCT DRAFT MODELS FROM KEGG 

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

%% 6. CONSTRUCT DRAFT MODEL FROM METACYC
%
modelAfmMC=getMetaCycModelForOrganism('afm','afm_prot.faa',false,false,false,100,55) 

%{

NOTE: Importing MetaCyc reactions from C:/Users/zorrilla/Documents/GitHub/RAVEN/external/metacyc/metaCycRxns.mat.
MetaCyc reactions loaded
NOTE: Importing MetaCyc enzymes and reaction-enzyme association from C:/Users/zorrilla/Documents/GitHub/RAVEN/external/metacyc/metaCycEnzymes.mat.
MetaCyc enzymes loaded
NOTE: Importing MetaCyc metabolites from C:/Users/zorrilla/Documents/GitHub/RAVEN/external/metacyc/metaCycMets.mat.
MetaCyc compounds loaded
The full MetaCyc model loaded
BLASTing "MetaCyc" against "afm"..
BLASTing "afm" against "MetaCyc"..
Completed searching against MetaCyc protein sequences.
NOTE: Importing MetaCyc metabolites from C:/Users/zorrilla/Documents/GitHub/RAVEN/external/metacyc/metaCycMets.mat.

modelAfmMC = 

  struct with fields:

                     id: 'afm'
            description: 'Generated by homology with MetaCyc database'
                   rxns: {1079×1 cell}
               rxnNames: {1079×1 cell}
                eccodes: {1079×1 cell}
             subSystems: {1079×1 cell}
             rxnMiriams: {1079×1 cell}
          rxnReferences: {1079×1 cell}
                     lb: [1079×1 double]
                     ub: [1079×1 double]
                    rev: [1079×1 double]
                      c: [1079×1 double]
              equations: {1079×1 cell}
                  genes: {1122×1 cell}
             rxnGeneMat: [1079×1122 double]
    rxnConfidenceScores: [1079×1 double]
                grRules: {1079×1 cell}
                      S: [1446×1079 double]
                   mets: {1446×1 cell}
               metNames: {1446×1 cell}
            metFormulas: {1446×1 cell}
             metCharges: [1446×1 double]
                      b: [1446×1 double]
                 inchis: {1446×1 cell}
             metMiriams: {1446×1 cell}
                  comps: {'s'}
              compNames: {'System'}
               metComps: [1446×1 double]

%}

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

%{

NOTE: A total of 1124 reactions in the KEGG model were mapped to MetaCyc.
NOTE: 490 reactions already in MetaCyc model, 634 will be combined.
NOTE: A shrinked KEGG model with 2122 reactions and 1738 metabolites was obtained.
NOTE: A total of 942 metabolites from the shrinked KEGG model were mapped to MetaCyc again.
NOTE: Importing KEGG metabolites from C:/Users/zorrilla/Documents/GitHub/RAVEN/external/kegg/keggMets.mat.

modelAfmMCKegg = 

  struct with fields:

                     id: 'COMBINED'
            description: 'Combined model from MetaCyc and KEGG draft models'
                   rxns: {3835×1 cell}
               rxnNames: {3835×1 cell}
                eccodes: {3835×1 cell}
             subSystems: {3835×1 cell}
             rxnMiriams: {3835×1 cell}
          rxnReferences: {3835×1 cell}
                     lb: [3835×1 double]
                     ub: [3835×1 double]
                    rev: [3835×1 double]
                      c: [3835×1 double]
              equations: {3835×1 cell}
                  genes: {3107×1 cell}
                grRules: {3835×1 cell}
                   mets: {2890×1 cell}
             metCharges: [2890×1 double]
                      b: [2890×1 double]
                 inchis: {2890×1 cell}
             metMiriams: {2890×1 cell}
                  comps: {'s'}
              compNames: {'System'}
               metComps: [2890×1 double]
                rxnFrom: {3835×1 cell}
               geneFrom: {3107×1 cell}
            grRulesKEGG: {3835×1 cell}
                      S: [2890×3835 double]
                metFrom: {2890×1 cell}
               metNames: {2890×1 cell}
            metFormulas: {2890×1 cell}
    rxnConfidenceScores: [3835×1 double]
             rxnGeneMat: [3835×3107 double]
%}

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

%{

modelAfmHomologyMergeComps = 

  struct with fields:

                     id: 'afm'
            description: 'Generated by getModelFromHomology using pch'
             annotation: [1×1 struct]
                   rxns: {1128×1 cell}
                   mets: {823×1 cell}
                      S: [823×1128 double]
                     lb: [1128×1 double]
                     ub: [1128×1 double]
                    rev: [1128×1 double]
                      c: [1128×1 double]
                      b: [823×1 double]
                  comps: {'s'}
              compNames: {'System'}
            compOutside: {''}
            compMiriams: {5×1 cell}
               rxnNames: {1128×1 cell}
                grRules: {1128×1 cell}
             rxnGeneMat: [1128×975 double]
             subSystems: {1128×1 cell}
                eccodes: {1128×1 cell}
                  genes: {975×1 cell}
              geneComps: [975×1 double]
               metNames: {823×1 cell}
               metComps: [823×1 double]
                 inchis: {823×1 cell}
            metFormulas: {823×1 cell}
             metMiriams: {823×1 cell}
               rxnNotes: {1128×1 cell}
    rxnConfidenceScores: [1128×1 double]

%}

%fill gaps before localization, merge comps of Pch model before fill gaps

%[modelPchMerged, deletedRxns, duplicateRxns]=mergeCompartments(modelPch,true);

%{

modelPchMerged = 

  struct with fields:

                id: 'pch'
       description: 'Penicillium chrysogenum genome-scale model'
        annotation: [1×1 struct]
              rxns: {1231×1 cell}
              mets: {848×1 cell}
                 S: [848×1231 double]
                lb: [1231×1 double]
                ub: [1231×1 double]
               rev: [1231×1 double]
                 c: [1231×1 double]
                 b: [848×1 double]
             comps: {'s'}
         compNames: {'System'}
       compOutside: {''}
       compMiriams: {5×1 cell}
          rxnNames: {1231×1 cell}
           grRules: {1231×1 cell}
        rxnGeneMat: [1231×948 double]
        subSystems: {1231×1 cell}
           eccodes: {1231×1 cell}
             genes: {948×1 cell}
         geneComps: [948×1 double]
       geneMiriams: {948×1 cell}
    geneShortNames: {948×1 cell}
          metNames: {848×1 cell}
          metComps: [848×1 double]
            inchis: {848×1 cell}
       metFormulas: {848×1 cell}
        metMiriams: {848×1 cell}

%}

[newConnected, cannotConnect, addedRxns, newModel, exitFlag]=fillGaps(modelAfmHomologyMergeComps,{modelPchMerged})

%{

MILP detected.
Academic license - for non-commercial use only
Optimize a model with 2264 rows, 6530 columns and 15686 nonzeros
Variable types: 5072 continuous, 1458 integer (0 binary)
Coefficient statistics:
  Matrix range     [1e-05, 3e+03]
  Objective range  [1e+00, 1e+00]
  Bounds range     [1e-01, 1e+03]
  RHS range        [0e+00, 0e+00]
Presolve removed 721 rows and 3654 columns
Presolve time: 0.09s
Presolved: 1543 rows, 2876 columns, 8659 nonzeros
Variable types: 1871 continuous, 1005 integer (1005 binary)
Presolve removed 1008 rows and 1018 columns
Presolved: 535 rows, 1858 columns, 6596 nonzeros


Root relaxation: objective 6.400060e+01, 395 iterations, 0.03 seconds

    Nodes    |    Current Node    |     Objective Bounds      |     Work
 Expl Unexpl |  Obj  Depth IntInf | Incumbent    BestBd   Gap | It/Node Time

     0     0   64.00060    0    4          -   64.00060      -     -    0s
H    0     0                      69.0000000   64.00060  7.25%     -    0s
H    0     0                      68.0000005   64.00060  5.88%     -    0s
     0     0     cutoff    0        68.00000   68.00000  0.00%     -    0s

Cutting planes:
  Cover: 2
  Implied bound: 2
  Clique: 1
  Flow cover: 1

Explored 1 nodes (623 simplex iterations) in 0.33 seconds
Thread count was 4 (of 4 available processors)

Solution count 2: 68 69 

Optimal solution found (tolerance 1.00e-09)
Best objective 6.800000051886e+01, best bound 6.800000051886e+01, gap 0.0000%

newConnected =

  141×1 cell array

    {'2ohpengOUT' }
    {'4hpoaIN'    }
    {'4ohpenvOUT' }
    {'6apaOUT'    }
    {'8hpaOUT'    }
    {'acvOUT'     }
    {'adIN'       }
    {'amiaceOUT'  }
    {'anIN'       }
    {'arabinIN'   }
    {'bmOUT'      }
    {'c181IN'     }
    {'c182IN'     }
    {'cyneIN'     }
    {'cytsIN'     }
    {'dglcIN'     }
    {'glcn15lacIN'}
    {'gnIN'       }
    {'hisIN'      }
    {'hyxnIN'     }
    {'ipnOUT'     }
    {'opcOUT'     }
    {'paaIN'      }
    {'penartOUT'  }
    {'pendhOUT'   }
    {'pengOUT'    }
    {'pengaOUT'   }
    {'penkOUT'    }
    {'penvOUT'    }
    {'penvaOUT'   }
    {'pheIN'      }
    {'piIN'       }
    {'poaIN'      }
    {'proIN'      }
    {'proteinOUT' }
    {'r0051'      }
    {'r0065'      }
    {'r0066'      }
    {'r0069'      }
    {'r0121'      }
    {'r0129'      }
    {'r0130'      }
    {'r0132'      }
    {'r0148'      }
    {'r0151'      }
    {'r0152'      }
    {'r0155'      }
    {'r0194'      }
    {'r0195'      }
    {'r0202'      }
    {'r0233'      }
    {'r0243'      }
    {'r0245'      }
    {'r0248'      }
    {'r0255'      }
    {'r0256'      }
    {'r0299'      }
    {'r0310'      }
    {'r0311'      }
    {'r0334'      }
    {'r0335'      }
    {'r0343'      }
    {'r0344'      }
    {'r0346'      }
    {'r0358'      }
    {'r0359'      }
    {'r0360'      }
    {'r0361'      }
    {'r0409'      }
    {'r0413'      }
    {'r0425'      }
    {'r0449'      }
    {'r0497'      }
    {'r0498'      }
    {'r0499'      }
    {'r0500'      }
    {'r0501'      }
    {'r0536'      }
    {'r0559'      }
    {'r0614'      }
    {'r0626'      }
    {'r0627'      }
    {'r0638'      }
    {'r0639'      }
    {'r0640'      }
    {'r0642'      }
    {'r0643'      }
    {'r0644'      }
    {'r0645'      }
    {'r0647'      }
    {'r0658'      }
    {'r0701'      }
    {'r0716'      }
    {'r0724'      }
    {'r0725'      }
    {'r0743'      }
    {'r0744'      }
    {'r0747'      }
    {'r0748'      }
    {'r0749'      }
    {'r0750'      }
    {'r0763'      }
    {'r0764'      }
    {'r0766'      }
    {'r0769'      }
    {'r0770'      }
    {'r0771'      }
    {'r0772'      }
    {'r0774'      }
    {'r0776'      }
    {'r0777'      }
    {'r0778'      }
    {'r0780'      }
    {'r0781'      }
    {'r0783'      }
    {'r0784'      }
    {'r0785'      }
    {'r0786'      }
    {'r0788'      }
    {'r0790'      }
    {'r0791'      }
    {'r0796'      }
    {'r0806'      }
    {'r0807'      }
    {'r0808'      }
    {'r0823'      }
    {'r0825'      }
    {'r1124'      }
    {'r1125'      }
    {'r1139'      }
    {'r1140'      }
    {'r1141'      }
    {'r1143'      }
    {'r1144'      }
    {'r1152'      }
    {'r1162'      }
    {'r1454'      }
    {'trpIN'      }
    {'uraIN'      }
    {'uriIN'      }
    {'xanIN'      }


cannotConnect =

  250×1 cell array

    {'balaIN' }
    {'c161IN' }
    {'c162IN' }
    {'c171IN' }
    {'c183IN' }
    {'c4odOUT'}
    {'fmnIN'  }
    {'myoiIN' }
    {'nicaIN' }
    {'nicdIN' }
    {'pabaIN' }
    {'pimIN'  }
    {'quinIN' }
    {'r0035'  }
    {'r0036'  }
    {'r0040'  }
    {'r0050'  }
    {'r0061'  }
    {'r0071'  }
    {'r0108'  }
    {'r0109'  }
    {'r0135'  }
    {'r0136'  }
    {'r0137'  }
    {'r0138'  }
    {'r0139'  }
    {'r0140'  }
    {'r0141'  }
    {'r0142'  }
    {'r0143'  }
    {'r0144'  }
    {'r0145'  }
    {'r0146'  }
    {'r0156'  }
    {'r0162'  }
    {'r0169'  }
    {'r0170'  }
    {'r0174'  }
    {'r0180'  }
    {'r0213'  }
    {'r0224'  }
    {'r0234'  }
    {'r0235'  }
    {'r0238'  }
    {'r0259'  }
    {'r0260'  }
    {'r0261'  }
    {'r0262'  }
    {'r0263'  }
    {'r0264'  }
    {'r0265'  }
    {'r0266'  }
    {'r0267'  }
    {'r0268'  }
    {'r0269'  }
    {'r0270'  }
    {'r0272'  }
    {'r0273'  }
    {'r0274'  }
    {'r0275'  }
    {'r0280'  }
    {'r0281'  }
    {'r0282'  }
    {'r0283'  }
    {'r0284'  }
    {'r0287'  }
    {'r0288'  }
    {'r0289'  }
    {'r0301'  }
    {'r0315'  }
    {'r0329'  }
    {'r0340'  }
    {'r0341'  }
    {'r0342'  }
    {'r0345'  }
    {'r0347'  }
    {'r0348'  }
    {'r0356'  }
    {'r0364'  }
    {'r0365'  }
    {'r0366'  }
    {'r0371'  }
    {'r0385'  }
    {'r0386'  }
    {'r0391'  }
    {'r0392'  }
    {'r0393'  }
    {'r0394'  }
    {'r0397'  }
    {'r0402'  }
    {'r0405'  }
    {'r0407'  }
    {'r0408'  }
    {'r0410'  }
    {'r0412'  }
    {'r0420'  }
    {'r0424'  }
    {'r0427'  }
    {'r0428'  }
    {'r0429'  }
    {'r0432'  }
    {'r0433'  }
    {'r0434'  }
    {'r0435'  }
    {'r0437'  }
    {'r0438'  }
    {'r0445'  }
    {'r0446'  }
    {'r0447'  }
    {'r0448'  }
    {'r0451'  }
    {'r0453'  }
    {'r0454'  }
    {'r0455'  }
    {'r0456'  }
    {'r0457'  }
    {'r0458'  }
    {'r0459'  }
    {'r0461'  }
    {'r0462'  }
    {'r0463'  }
    {'r0464'  }
    {'r0467'  }
    {'r0468'  }
    {'r0469'  }
    {'r0471'  }
    {'r0472'  }
    {'r0473'  }
    {'r0474'  }
    {'r0475'  }
    {'r0476'  }
    {'r0477'  }
    {'r0478'  }
    {'r0479'  }
    {'r0480'  }
    {'r0481'  }
    {'r0482'  }
    {'r0483'  }
    {'r0486'  }
    {'r0487'  }
    {'r0488'  }
    {'r0491'  }
    {'r0496'  }
    {'r0502'  }
    {'r0503'  }
    {'r0506'  }
    {'r0510'  }
    {'r0511'  }
    {'r0512'  }
    {'r0514'  }
    {'r0515'  }
    {'r0516'  }
    {'r0517'  }
    {'r0521'  }
    {'r0523'  }
    {'r0525'  }
    {'r0526'  }
    {'r0527'  }
    {'r0529'  }
    {'r0531'  }
    {'r0532'  }
    {'r0533'  }
    {'r0534'  }
    {'r0535'  }
    {'r0537'  }
    {'r0538'  }
    {'r0539'  }
    {'r0561'  }
    {'r0564'  }
    {'r0565'  }
    {'r0573'  }
    {'r0576'  }
    {'r0585'  }
    {'r0586'  }
    {'r0597'  }
    {'r0602'  }
    {'r0603'  }
    {'r0604'  }
    {'r0607'  }
    {'r0608'  }
    {'r0610'  }
    {'r0611'  }
    {'r0612'  }
    {'r0630'  }
    {'r0631'  }
    {'r0641'  }
    {'r0659'  }
    {'r0660'  }
    {'r0661'  }
    {'r0691'  }
    {'r0697'  }
    {'r0707'  }
    {'r0710'  }
    {'r0711'  }
    {'r0712'  }
    {'r0717'  }
    {'r0723'  }
    {'r0726'  }
    {'r0728'  }
    {'r0729'  }
    {'r0730'  }
    {'r0731'  }
    {'r0746'  }
    {'r0752'  }
    {'r0753'  }
    {'r0755'  }
    {'r0757'  }
    {'r0759'  }
    {'r0762'  }
    {'r0824'  }
    {'r0826'  }
    {'r0827'  }
    {'r0828'  }
    {'r1029'  }
    {'r1032'  }
    {'r1035'  }
    {'r1041'  }
    {'r1046'  }
    {'r1051'  }
    {'r1056'  }
    {'r1061'  }
    {'r1067'  }
    {'r1072'  }
    {'r1077'  }
    {'r1082'  }
    {'r1087'  }
    {'r1096'  }
    {'r1097'  }
    {'r1104'  }
    {'r1105'  }
    {'r1118'  }
    {'r1119'  }
    {'r1120'  }
    {'r1121'  }
    {'r1122'  }
    {'r1123'  }
    {'r1126'  }
    {'r1127'  }
    {'r1128'  }
    {'r1129'  }
    {'r1131'  }
    {'r1132'  }
    {'r1133'  }
    {'r1134'  }
    {'r1135'  }
    {'r1138'  }
    {'r1142'  }
    {'r1145'  }
    {'r1146'  }
    {'thmIN'  }


addedRxns =

  69×1 cell array

    {'r0001'}
    {'r0044'}
    {'r0118'}
    {'r0154'}
    {'r0193'}
    {'r0236'}
    {'r0237'}
    {'r0244'}
    {'r0249'}
    {'r0250'}
    {'r0253'}
    {'r0258'}
    {'r0271'}
    {'r0276'}
    {'r0277'}
    {'r0278'}
    {'r0279'}
    {'r0285'}
    {'r0286'}
    {'r0300'}
    {'r0328'}
    {'r0404'}
    {'r0406'}
    {'r0466'}
    {'r0489'}
    {'r0490'}
    {'r0492'}
    {'r0504'}
    {'r0505'}
    {'r0507'}
    {'r0508'}
    {'r0509'}
    {'r0513'}
    {'r0609'}
    {'r0615'}
    {'r0637'}
    {'r0646'}
    {'r0649'}
    {'r0745'}
    {'r0800'}
    {'r0801'}
    {'r0802'}
    {'r0803'}
    {'r0809'}
    {'r0810'}
    {'r0812'}
    {'r0813'}
    {'r0814'}
    {'r0816'}
    {'r0817'}
    {'r0818'}
    {'r0819'}
    {'r0820'}
    {'r1153'}
    {'r1163'}
    {'r1453'}
    {'r1455'}
    {'r1456'}
    {'r1457'}
    {'r1458'}
    {'r1459'}
    {'r1460'}
    {'r1461'}
    {'r1462'}
    {'r1463'}
    {'r1464'}
    {'r1465'}
    {'r1466'}
    {'r1467'}


newModel = 

  struct with fields:

                     id: 'MERGED'
            description: ''
             annotation: [1×1 struct]
                   rxns: {1197×1 cell}
                   mets: {843×1 cell}
                      S: [843×1197 double]
                     lb: [1197×1 double]
                     ub: [1197×1 double]
                    rev: [1197×1 double]
                      c: [1197×1 double]
                      b: [843×1 double]
                  comps: {'s'}
              compNames: {'System'}
            compOutside: {''}
            compMiriams: {5×1 cell}
               rxnNames: {1197×1 cell}
                grRules: {1197×1 cell}
             rxnGeneMat: [1197×1923 double]
             subSystems: {1197×1 cell}
                eccodes: {1197×1 cell}
                  genes: {1923×1 cell}
              geneComps: [1923×1 double]
               metNames: {843×1 cell}
               metComps: [843×1 double]
                 inchis: {843×1 cell}
            metFormulas: {843×1 cell}
             metMiriams: {843×1 cell}
               rxnNotes: {1197×1 cell}
    rxnConfidenceScores: [1197×1 double]
                rxnFrom: {1197×1 cell}
                metFrom: {843×1 cell}
               geneFrom: {1923×1 cell}
         geneShortNames: {1923×1 cell}
            geneMiriams: {1923×1 cell}


exitFlag =

     1

%}

[modelAfmHomologyDeeploc, geneLocalization, transportStruct, score, removedRxns] = predictLocalization(modelAfmHomology,GSSsubset,'Cytoplasm')

%{

WARNING: The model structure contains information about gene compartmentalization. This is not supported by this function. The geneComps field has been deleted

WARNING: Reaction r0088 contains nested and/or-relations. Large risk of errors

WARNING: Reaction r0093 contains nested and/or-relations. Large risk of errors

WARNING: Reaction r0094 contains nested and/or-relations. Large risk of errors

WARNING: Reaction r0095 contains nested and/or-relations. Large risk of errors

WARNING: Reaction r0113 contains nested and/or-relations. Large risk of errors

WARNING: Reaction r0122 contains nested and/or-relations. Large risk of errors

WARNING: Reaction r0409 contains nested and/or-relations. Large risk of errors

WARNING: Reaction r0425 contains nested and/or-relations. Large risk of errors

WARNING: Reaction r0465 contains nested and/or-relations. Large risk of errors

WARNING: Reaction r0631 contains nested and/or-relations. Large risk of errors

WARNING: Reaction r0651 contains nested and/or-relations. Large risk of errors

WARNING: Reaction r0673 contains nested and/or-relations. Large risk of errors

WARNING: Reaction r0677 contains nested and/or-relations. Large risk of errors

WARNING: Reaction r0702 contains nested and/or-relations. Large risk of errors

WARNING: Reaction r0705 contains nested and/or-relations. Large risk of errors


modelAfmHomologyDeeploc = 

  struct with fields:

                     id: 'afm'
            description: 'Generated by getModelFromHomology using pch'
             annotation: [1×1 struct]
                   rxns: {1966×1 cell}
                   mets: {1035×1 cell}
                      S: [1035×1966 double]
                     lb: [1966×1 double]
                     ub: [1966×1 double]
                    rev: [1966×1 double]
                      c: [1966×1 double]
                      b: [1035×1 double]
                  comps: {4×1 cell}
              compNames: {4×1 cell}
            compOutside: {4×1 cell}
            compMiriams: {5×1 cell}
               rxnNames: {1966×1 cell}
                grRules: {1966×1 cell}
             rxnGeneMat: [1966×800 double]
             subSystems: {1966×1 cell}
                eccodes: {1966×1 cell}
                  genes: {800×1 cell}
               metNames: {1035×1 cell}
               metComps: [1035×1 double]
                 inchis: {1035×1 cell}
            metFormulas: {1035×1 cell}
             metMiriams: {1035×1 cell}
               rxnNotes: {1966×1 cell}
    rxnConfidenceScores: [1966×1 double]


geneLocalization = 

  struct with fields:

    genes: {800×1 cell}
    comps: {800×1 cell}


transportStruct = 

  struct with fields:

      mets: {'L-tryptophan'}
    toComp: {'Mitochondrion'}


score = 

  struct with fields:

     totScore: 398
    geneScore: 400
    transCost: 2


removedRxns =

  352×1 cell array

    {'r0065'       }
    {'r0066'       }
    {'r0071'       }
    {'r0121'       }
    {'r0155'       }
    {'r0162'       }
    {'r0162'       }
    {'r0167'       }
    {'r0168'       }
    {'r0169'       }
    {'r0170'       }
    {'r0172'       }
    {'r0174'       }
    {'r0176'       }
    {'r0194'       }
    {'r0200'       }
    {'r0202'       }
    {'r0207'       }
    {'r0213'       }
    {'r0220'       }
    {'r0239'       }
    {'r0255'       }
    {'r0256'       }
    {'r0284'       }
    {'r0284'       }
    {'r0287'       }
    {'r0301'       }
    {'r0310'       }
    {'r0311'       }
    {'r0318'       }
    {'r0320'       }
    {'r0326'       }
    {'r0334'       }
    {'r0335'       }
    {'r0356'       }
    {'r0405'       }
    {'r0449'       }
    {'r0457'       }
    {'r0463'       }
    {'r0487'       }
    {'r0506'       }
    {'r0510'       }
    {'r0514'       }
    {'r0515'       }
    {'r0516'       }
    {'r0523'       }
    {'r0534'       }
    {'r0537'       }
    {'r0539'       }
    {'r0573'       }
    {'r0597'       }
    {'r0602'       }
    {'r0603'       }
    {'r0604'       }
    {'r0610'       }
    {'r0627'       }
    {'r0630'       }
    {'r0638'       }
    {'r0726'       }
    {'r0728'       }
    {'r0743'       }
    {'r0749'       }
    {'r0752'       }
    {'r0759'       }
    {'r0763'       }
    {'r0764'       }
    {'r0766'       }
    {'r0769'       }
    {'r0770'       }
    {'r0771'       }
    {'r0772'       }
    {'r0774'       }
    {'r0776'       }
    {'r0777'       }
    {'r0778'       }
    {'r0780'       }
    {'r0781'       }
    {'r0783'       }
    {'r0784'       }
    {'r0785'       }
    {'r0786'       }
    {'r0788'       }
    {'r0790'       }
    {'r0791'       }
    {'r0796'       }
    {'r0806'       }
    {'r1152'       }
    {'r1158'       }
    {'r1159'       }
    {'r1161'       }
    {'r1162'       }
    {'r0192'       }
    {'r0241'       }
    {'r0065_EXP_2' }
    {'r0065_EXP_3' }
    {'r0066_EXP_2' }
    {'r0066_EXP_3' }
    {'r0066_EXP_4' }
    {'r0071_EXP_2' }
    {'r0121_EXP_2' }
    {'r0167_EXP_2' }
    {'r0167_EXP_3' }
    {'r0167_EXP_4' }
    {'r0167_EXP_5' }
    {'r0167_EXP_6' }
    {'r0167_EXP_7' }
    {'r0167_EXP_8' }
    {'r0167_EXP_9' }
    {'r0167_EXP_10'}
    {'r0167_EXP_11'}
    {'r0167_EXP_12'}
    {'r0168_EXP_2' }
    {'r0168_EXP_3' }
    {'r0168_EXP_4' }
    {'r0168_EXP_5' }
    {'r0168_EXP_6' }
    {'r0168_EXP_7' }
    {'r0168_EXP_8' }
    {'r0168_EXP_9' }
    {'r0168_EXP_10'}
    {'r0168_EXP_11'}
    {'r0168_EXP_12'}
    {'r0169_EXP_2' }
    {'r0169_EXP_3' }
    {'r0169_EXP_4' }
    {'r0169_EXP_5' }
    {'r0169_EXP_6' }
    {'r0169_EXP_7' }
    {'r0169_EXP_8' }
    {'r0169_EXP_9' }
    {'r0169_EXP_10'}
    {'r0170_EXP_2' }
    {'r0170_EXP_3' }
    {'r0170_EXP_4' }
    {'r0170_EXP_5' }
    {'r0170_EXP_6' }
    {'r0170_EXP_7' }
    {'r0170_EXP_8' }
    {'r0170_EXP_9' }
    {'r0170_EXP_10'}
    {'r0174_EXP_2' }
    {'r0176_EXP_2' }
    {'r0176_EXP_3' }
    {'r0176_EXP_4' }
    {'r0194_EXP_2' }
    {'r0194_EXP_3' }
    {'r0194_EXP_4' }
    {'r0194_EXP_5' }
    {'r0194_EXP_6' }
    {'r0200_EXP_2' }
    {'r0200_EXP_3' }
    {'r0200_EXP_4' }
    {'r0200_EXP_5' }
    {'r0200_EXP_6' }
    {'r0200_EXP_7' }
    {'r0200_EXP_8' }
    {'r0200_EXP_9' }
    {'r0200_EXP_10'}
    {'r0202_EXP_2' }
    {'r0202_EXP_3' }
    {'r0202_EXP_4' }
    {'r0202_EXP_5' }
    {'r0207_EXP_2' }
    {'r0207_EXP_3' }
    {'r0207_EXP_4' }
    {'r0207_EXP_5' }
    {'r0207_EXP_6' }
    {'r0207_EXP_7' }
    {'r0207_EXP_8' }
    {'r0207_EXP_9' }
    {'r0213_EXP_2' }
    {'r0213_EXP_3' }
    {'r0213_EXP_4' }
    {'r0220_EXP_2' }
    {'r0220_EXP_3' }
    {'r0220_EXP_4' }
    {'r0220_EXP_5' }
    {'r0220_EXP_6' }
    {'r0220_EXP_7' }
    {'r0220_EXP_8' }
    {'r0220_EXP_9' }
    {'r0239_EXP_2' }
    {'r0239_EXP_3' }
    {'r0239_EXP_4' }
    {'r0239_EXP_5' }
    {'r0255_EXP_2' }
    {'r0255_EXP_3' }
    {'r0255_EXP_4' }
    {'r0255_EXP_5' }
    {'r0256_EXP_2' }
    {'r0256_EXP_3' }
    {'r0310_EXP_2' }
    {'r0310_EXP_3' }
    {'r0311_EXP_2' }
    {'r0311_EXP_3' }
    {'r0311_EXP_4' }
    {'r0318_EXP_2' }
    {'r0318_EXP_3' }
    {'r0318_EXP_4' }
    {'r0320_EXP_2' }
    {'r0320_EXP_3' }
    {'r0320_EXP_4' }
    {'r0487_EXP_2' }
    {'r0506_EXP_2' }
    {'r0523_EXP_2' }
    {'r0534_EXP_2' }
    {'r0534_EXP_3' }
    {'r0573_EXP_2' }
    {'r0597_EXP_2' }
    {'r0627_EXP_2' }
    {'r0627_EXP_3' }
    {'r0627_EXP_4' }
    {'r0627_EXP_5' }
    {'r0627_EXP_6' }
    {'r0627_EXP_7' }
    {'r0638_EXP_2' }
    {'r0726_EXP_2' }
    {'r0726_EXP_3' }
    {'r0726_EXP_4' }
    {'r0726_EXP_5' }
    {'r0726_EXP_6' }
    {'r0726_EXP_7' }
    {'r0726_EXP_8' }
    {'r0726_EXP_9' }
    {'r0726_EXP_10'}
    {'r0726_EXP_11'}
    {'r0728_EXP_2' }
    {'r0743_EXP_2' }
    {'r0743_EXP_3' }
    {'r0749_EXP_2' }
    {'r0752_EXP_2' }
    {'r0752_EXP_3' }
    {'r0759_EXP_2' }
    {'r0764_EXP_2' }
    {'r0764_EXP_3' }
    {'r0764_EXP_4' }
    {'r0766_EXP_2' }
    {'r0770_EXP_2' }
    {'r0772_EXP_2' }
    {'r0774_EXP_2' }
    {'r0776_EXP_2' }
    {'r0777_EXP_2' }
    {'r0778_EXP_2' }
    {'r0781_EXP_2' }
    {'r0784_EXP_2' }
    {'r0785_EXP_2' }
    {'r0788_EXP_2' }
    {'r0790_EXP_2' }
    {'r0806_EXP_2' }
    {'r0806_EXP_3' }
    {'r0806_EXP_4' }
    {'r1161_EXP_2' }
    {'r1161_EXP_3' }
    {'r1161_EXP_4' }
    {'r1161_EXP_5' }
    {'r0192_EXP_2' }
    {'r0192_EXP_3' }
    {'r0192_EXP_4' }
    {'r0241_EXP_2' }
    {'r0241_EXP_3' }
    {'r0241_EXP_4' }
    {'r0050'       }
    {'r0050'       }
    {'r0051'       }
    {'r0069'       }
    {'r0148'       }
    {'r0159'       }
    {'r0160'       }
    {'r0177'       }
    {'r0231'       }
    {'r0288'       }
    {'r0316'       }
    {'r0332'       }
    {'r0468'       }
    {'r0488'       }
    {'r0511'       }
    {'r0517'       }
    {'r0517'       }
    {'r0536'       }
    {'r0538'       }
    {'r0611'       }
    {'r0631'       }
    {'r0639'       }
    {'r0729'       }
    {'r0744'       }
    {'r0807'       }
    {'r1157'       }
    {'r0051_EXP_2' }
    {'r0069_EXP_2' }
    {'r0159_EXP_2' }
    {'r0159_EXP_3' }
    {'r0159_EXP_4' }
    {'r0159_EXP_5' }
    {'r0159_EXP_6' }
    {'r0159_EXP_7' }
    {'r0160_EXP_2' }
    {'r0160_EXP_3' }
    {'r0160_EXP_4' }
    {'r0160_EXP_5' }
    {'r0160_EXP_6' }
    {'r0160_EXP_7' }
    {'r0177_EXP_2' }
    {'r0177_EXP_3' }
    {'r0177_EXP_4' }
    {'r0231_EXP_2' }
    {'r0316_EXP_2' }
    {'r0332_EXP_2' }
    {'r0488_EXP_2' }
    {'r0536_EXP_2' }
    {'r0611_EXP_2' }
    {'r0631_EXP_2' }
    {'r0639_EXP_2' }
    {'r0729_EXP_2' }
    {'r0744_EXP_2' }
    {'r0744_EXP_3' }
    {'r0744_EXP_4' }
    {'r0744_EXP_5' }
    {'r0744_EXP_6' }
    {'r0744_EXP_7' }
    {'r0807_EXP_2' }
    {'r0807_EXP_3' }
    {'r0807_EXP_4' }
    {'r0188'       }
    {'r0289'       }
    {'r0330'       }
    {'r0333'       }
    {'r0469'       }
    {'r0512'       }
    {'r0642'       }
    {'r0188_EXP_2' }
    {'r0188_EXP_3' }
    {'r0188_EXP_4' }
    {'r0188_EXP_5' }
    {'r0188_EXP_6' }
    {'r0188_EXP_7' }
    {'r0324'       }
    {'r0467'       }
    {'r0640'       }
    {'r0324_EXP_2' }
    {'r0324_EXP_3' }
    {'r0324_EXP_4' }
    {'r0324_EXP_5' }
    {'r0324_EXP_6' }
    {'r0464'       }
    {'r0643'       }
    {'r0644'       }
    {'r0645'       }
    {'r0647'       }
    {'r0647_EXP_2' }
    {'r0641'       }
    {'r0641_EXP_2' }
    {'r0641_EXP_3' }


%}

modelAfmHomologyDeeploc = rmfield(modelAfmHomologyDeeploc,'geneComps');
exportToExcelFormat(modelAfmHomologyDeeploc,'modelAfmHomologyDeeploc.xlsx'); 

% try contracting model, gap
% fill using unmergedPch (careful with naming)

% try same thing but with the modelAfmHomologyPreMerged

%[modelAfmHomologyPreMergedCello, geneLocalization, transportStruct, score, removedRxns] = predictLocalization(modelAfmHomologyPreMerged ,GSSsubset,'cytop')

%{

WARNING: The model structure contains information about gene compartmentalization. This is not supported by this function. The geneComps field has been deleted

WARNING: Reaction r0088 contains nested and/or-relations. Large risk of errors

WARNING: Reaction r0093 contains nested and/or-relations. Large risk of errors

WARNING: Reaction r0094 contains nested and/or-relations. Large risk of errors

WARNING: Reaction r0095 contains nested and/or-relations. Large risk of errors

WARNING: Reaction r0113 contains nested and/or-relations. Large risk of errors

WARNING: Reaction r0122 contains nested and/or-relations. Large risk of errors

WARNING: Reaction r0409 contains nested and/or-relations. Large risk of errors

WARNING: Reaction r0425 contains nested and/or-relations. Large risk of errors

WARNING: Reaction r0465 contains nested and/or-relations. Large risk of errors

WARNING: Reaction r0631 contains nested and/or-relations. Large risk of errors

WARNING: Reaction r0651 contains nested and/or-relations. Large risk of errors

WARNING: Reaction r0673 contains nested and/or-relations. Large risk of errors

WARNING: Reaction r0677 contains nested and/or-relations. Large risk of errors

WARNING: Reaction r0702 contains nested and/or-relations. Large risk of errors

WARNING: Reaction r0705 contains nested and/or-relations. Large risk of errors


modelAfmHomologyPreMergedCello = 

  struct with fields:

                     id: 'afm'
            description: 'Generated by getModelFromHomology using pch'
             annotation: [1×1 struct]
                   rxns: {1966×1 cell}
                   mets: {906×1 cell}
                      S: [906×1966 double]
                     lb: [1966×1 double]
                     ub: [1966×1 double]
                    rev: [1966×1 double]
                      c: [1966×1 double]
                      b: [906×1 double]
                  comps: {4×1 cell}
              compNames: {4×1 cell}
            compOutside: {4×1 cell}
            compMiriams: {5×1 cell}
               rxnNames: {1966×1 cell}
                grRules: {1966×1 cell}
             rxnGeneMat: [1966×800 double]
             subSystems: {1966×1 cell}
                eccodes: {1966×1 cell}
                  genes: {800×1 cell}
               metNames: {906×1 cell}
               metComps: [906×1 double]
                 inchis: {906×1 cell}
            metFormulas: {906×1 cell}
             metMiriams: {906×1 cell}
               rxnNotes: {1966×1 cell}
    rxnConfidenceScores: [1966×1 double]


geneLocalization = 

  struct with fields:

    genes: {800×1 cell}
    comps: {800×1 cell}


transportStruct = 

  struct with fields:

      mets: {'N2-acetyl-L-ornithine'}
    toComp: {'mito'}


score = 

  struct with fields:

     totScore: 1.3988
    geneScore: 548.5844
    transCost: 2


removedRxns =

  352×1 cell array

    {'r0065'       }
    {'r0066'       }
    {'r0071'       }
    {'r0121'       }
    {'r0155'       }
    {'r0162'       }
    {'r0162'       }
    {'r0167'       }
    {'r0168'       }
    {'r0169'       }
    {'r0170'       }
    {'r0172'       }
    {'r0174'       }
    {'r0176'       }
    {'r0194'       }
    {'r0200'       }
    {'r0202'       }
    {'r0207'       }
    {'r0213'       }
    {'r0220'       }
    {'r0239'       }
    {'r0255'       }
    {'r0256'       }
    {'r0284'       }
    {'r0284'       }
    {'r0287'       }
    {'r0301'       }
    {'r0310'       }
    {'r0311'       }
    {'r0318'       }
    {'r0320'       }
    {'r0326'       }
    {'r0334'       }
    {'r0335'       }
    {'r0356'       }
    {'r0405'       }
    {'r0449'       }
    {'r0457'       }
    {'r0463'       }
    {'r0487'       }
    {'r0506'       }
    {'r0510'       }
    {'r0514'       }
    {'r0515'       }
    {'r0516'       }
    {'r0523'       }
    {'r0534'       }
    {'r0537'       }
    {'r0539'       }
    {'r0573'       }
    {'r0597'       }
    {'r0602'       }
    {'r0603'       }
    {'r0604'       }
    {'r0610'       }
    {'r0627'       }
    {'r0630'       }
    {'r0638'       }
    {'r0726'       }
    {'r0728'       }
    {'r0743'       }
    {'r0749'       }
    {'r0752'       }
    {'r0759'       }
    {'r0763'       }
    {'r0764'       }
    {'r0766'       }
    {'r0769'       }
    {'r0770'       }
    {'r0771'       }
    {'r0772'       }
    {'r0774'       }
    {'r0776'       }
    {'r0777'       }
    {'r0778'       }
    {'r0780'       }
    {'r0781'       }
    {'r0783'       }
    {'r0784'       }
    {'r0785'       }
    {'r0786'       }
    {'r0788'       }
    {'r0790'       }
    {'r0791'       }
    {'r0796'       }
    {'r0806'       }
    {'r1152'       }
    {'r1158'       }
    {'r1159'       }
    {'r1161'       }
    {'r1162'       }
    {'r0192'       }
    {'r0241'       }
    {'r0065_EXP_2' }
    {'r0065_EXP_3' }
    {'r0066_EXP_2' }
    {'r0066_EXP_3' }
    {'r0066_EXP_4' }
    {'r0071_EXP_2' }
    {'r0121_EXP_2' }
    {'r0167_EXP_2' }
    {'r0167_EXP_3' }
    {'r0167_EXP_4' }
    {'r0167_EXP_5' }
    {'r0167_EXP_6' }
    {'r0167_EXP_7' }
    {'r0167_EXP_8' }
    {'r0167_EXP_9' }
    {'r0167_EXP_10'}
    {'r0167_EXP_11'}
    {'r0167_EXP_12'}
    {'r0168_EXP_2' }
    {'r0168_EXP_3' }
    {'r0168_EXP_4' }
    {'r0168_EXP_5' }
    {'r0168_EXP_6' }
    {'r0168_EXP_7' }
    {'r0168_EXP_8' }
    {'r0168_EXP_9' }
    {'r0168_EXP_10'}
    {'r0168_EXP_11'}
    {'r0168_EXP_12'}
    {'r0169_EXP_2' }
    {'r0169_EXP_3' }
    {'r0169_EXP_4' }
    {'r0169_EXP_5' }
    {'r0169_EXP_6' }
    {'r0169_EXP_7' }
    {'r0169_EXP_8' }
    {'r0169_EXP_9' }
    {'r0169_EXP_10'}
    {'r0170_EXP_2' }
    {'r0170_EXP_3' }
    {'r0170_EXP_4' }
    {'r0170_EXP_5' }
    {'r0170_EXP_6' }
    {'r0170_EXP_7' }
    {'r0170_EXP_8' }
    {'r0170_EXP_9' }
    {'r0170_EXP_10'}
    {'r0174_EXP_2' }
    {'r0176_EXP_2' }
    {'r0176_EXP_3' }
    {'r0176_EXP_4' }
    {'r0194_EXP_2' }
    {'r0194_EXP_3' }
    {'r0194_EXP_4' }
    {'r0194_EXP_5' }
    {'r0194_EXP_6' }
    {'r0200_EXP_2' }
    {'r0200_EXP_3' }
    {'r0200_EXP_4' }
    {'r0200_EXP_5' }
    {'r0200_EXP_6' }
    {'r0200_EXP_7' }
    {'r0200_EXP_8' }
    {'r0200_EXP_9' }
    {'r0200_EXP_10'}
    {'r0202_EXP_2' }
    {'r0202_EXP_3' }
    {'r0202_EXP_4' }
    {'r0202_EXP_5' }
    {'r0207_EXP_2' }
    {'r0207_EXP_3' }
    {'r0207_EXP_4' }
    {'r0207_EXP_5' }
    {'r0207_EXP_6' }
    {'r0207_EXP_7' }
    {'r0207_EXP_8' }
    {'r0207_EXP_9' }
    {'r0213_EXP_2' }
    {'r0213_EXP_3' }
    {'r0213_EXP_4' }
    {'r0220_EXP_2' }
    {'r0220_EXP_3' }
    {'r0220_EXP_4' }
    {'r0220_EXP_5' }
    {'r0220_EXP_6' }
    {'r0220_EXP_7' }
    {'r0220_EXP_8' }
    {'r0220_EXP_9' }
    {'r0239_EXP_2' }
    {'r0239_EXP_3' }
    {'r0239_EXP_4' }
    {'r0239_EXP_5' }
    {'r0255_EXP_2' }
    {'r0255_EXP_3' }
    {'r0255_EXP_4' }
    {'r0255_EXP_5' }
    {'r0256_EXP_2' }
    {'r0256_EXP_3' }
    {'r0310_EXP_2' }
    {'r0310_EXP_3' }
    {'r0311_EXP_2' }
    {'r0311_EXP_3' }
    {'r0311_EXP_4' }
    {'r0318_EXP_2' }
    {'r0318_EXP_3' }
    {'r0318_EXP_4' }
    {'r0320_EXP_2' }
    {'r0320_EXP_3' }
    {'r0320_EXP_4' }
    {'r0487_EXP_2' }
    {'r0506_EXP_2' }
    {'r0523_EXP_2' }
    {'r0534_EXP_2' }
    {'r0534_EXP_3' }
    {'r0573_EXP_2' }
    {'r0597_EXP_2' }
    {'r0627_EXP_2' }
    {'r0627_EXP_3' }
    {'r0627_EXP_4' }
    {'r0627_EXP_5' }
    {'r0627_EXP_6' }
    {'r0627_EXP_7' }
    {'r0638_EXP_2' }
    {'r0726_EXP_2' }
    {'r0726_EXP_3' }
    {'r0726_EXP_4' }
    {'r0726_EXP_5' }
    {'r0726_EXP_6' }
    {'r0726_EXP_7' }
    {'r0726_EXP_8' }
    {'r0726_EXP_9' }
    {'r0726_EXP_10'}
    {'r0726_EXP_11'}
    {'r0728_EXP_2' }
    {'r0743_EXP_2' }
    {'r0743_EXP_3' }
    {'r0749_EXP_2' }
    {'r0752_EXP_2' }
    {'r0752_EXP_3' }
    {'r0759_EXP_2' }
    {'r0764_EXP_2' }
    {'r0764_EXP_3' }
    {'r0764_EXP_4' }
    {'r0766_EXP_2' }
    {'r0770_EXP_2' }
    {'r0772_EXP_2' }
    {'r0774_EXP_2' }
    {'r0776_EXP_2' }
    {'r0777_EXP_2' }
    {'r0778_EXP_2' }
    {'r0781_EXP_2' }
    {'r0784_EXP_2' }
    {'r0785_EXP_2' }
    {'r0788_EXP_2' }
    {'r0790_EXP_2' }
    {'r0806_EXP_2' }
    {'r0806_EXP_3' }
    {'r0806_EXP_4' }
    {'r1161_EXP_2' }
    {'r1161_EXP_3' }
    {'r1161_EXP_4' }
    {'r1161_EXP_5' }
    {'r0192_EXP_2' }
    {'r0192_EXP_3' }
    {'r0192_EXP_4' }
    {'r0241_EXP_2' }
    {'r0241_EXP_3' }
    {'r0241_EXP_4' }
    {'r0050'       }
    {'r0050'       }
    {'r0051'       }
    {'r0069'       }
    {'r0148'       }
    {'r0159'       }
    {'r0160'       }
    {'r0177'       }
    {'r0231'       }
    {'r0288'       }
    {'r0316'       }
    {'r0332'       }
    {'r0468'       }
    {'r0488'       }
    {'r0511'       }
    {'r0517'       }
    {'r0517'       }
    {'r0536'       }
    {'r0538'       }
    {'r0611'       }
    {'r0631'       }
    {'r0639'       }
    {'r0729'       }
    {'r0744'       }
    {'r0807'       }
    {'r1157'       }
    {'r0051_EXP_2' }
    {'r0069_EXP_2' }
    {'r0159_EXP_2' }
    {'r0159_EXP_3' }
    {'r0159_EXP_4' }
    {'r0159_EXP_5' }
    {'r0159_EXP_6' }
    {'r0159_EXP_7' }
    {'r0160_EXP_2' }
    {'r0160_EXP_3' }
    {'r0160_EXP_4' }
    {'r0160_EXP_5' }
    {'r0160_EXP_6' }
    {'r0160_EXP_7' }
    {'r0177_EXP_2' }
    {'r0177_EXP_3' }
    {'r0177_EXP_4' }
    {'r0231_EXP_2' }
    {'r0316_EXP_2' }
    {'r0332_EXP_2' }
    {'r0488_EXP_2' }
    {'r0536_EXP_2' }
    {'r0611_EXP_2' }
    {'r0631_EXP_2' }
    {'r0639_EXP_2' }
    {'r0729_EXP_2' }
    {'r0744_EXP_2' }
    {'r0744_EXP_3' }
    {'r0744_EXP_4' }
    {'r0744_EXP_5' }
    {'r0744_EXP_6' }
    {'r0744_EXP_7' }
    {'r0807_EXP_2' }
    {'r0807_EXP_3' }
    {'r0807_EXP_4' }
    {'r0188'       }
    {'r0289'       }
    {'r0330'       }
    {'r0333'       }
    {'r0469'       }
    {'r0512'       }
    {'r0642'       }
    {'r0188_EXP_2' }
    {'r0188_EXP_3' }
    {'r0188_EXP_4' }
    {'r0188_EXP_5' }
    {'r0188_EXP_6' }
    {'r0188_EXP_7' }
    {'r0324'       }
    {'r0467'       }
    {'r0640'       }
    {'r0324_EXP_2' }
    {'r0324_EXP_3' }
    {'r0324_EXP_4' }
    {'r0324_EXP_5' }
    {'r0324_EXP_6' }
    {'r0464'       }
    {'r0643'       }
    {'r0644'       }
    {'r0645'       }
    {'r0647'       }
    {'r0647_EXP_2' }
    {'r0641'       }
    {'r0641_EXP_2' }
    {'r0641_EXP_3' }

%}


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

