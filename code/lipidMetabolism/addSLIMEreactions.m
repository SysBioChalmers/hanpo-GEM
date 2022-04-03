function model = addSLIMEreactions(template, model, sourceModel)
% addSLIMEreactionss
%  Uses template reactions in a defined format to add SLIME reactions with
%  specified acyl-chain.
%
%  Input:
%    template    structure with defined fields:
%     .metName   strings with metabolite names, with placeholders to be
%                replaced with specified acyl chains.
%     .bbID      strings with metabolite id of relevant lipid backbone.
%     .bbMW      vector with molecular weights of lipid backbones without
%                acyl-chain.
%     .comps     strings with compartment abbreviation.
%     .chains    cell array with sets of acyl-chains, separated within
%                each set by comma, while individual sets of acyl-chains
%                are in individual cells.
%    model       model structure where reactions will be added to
%    sourcemodel model that will be queried to see if reaction already
%                exists, and in that case using the designated reaction id.
%
%  Output:
%    model       model structure with new lipid reactions added
%
%  Example of the template files can be found in the ComplementaryData/
%  reconstruction folder in github.com/sysbiochalmers/hanpo-GEM.
%  This function functions with models derived from yeast-GEM. Modify
%  acyl-chain and compartment information if used with different models.
%
%  Eduard Kerkhoven, 2019-07-03

dat.chain = {'16:0';      '16:1';         '18:0';     '18:1';   '18:2';      '18:3'};
dat.name  = {'palmitate'; 'palmitoleate'; 'stearate'; 'oleate'; 'linoleate'; 'linolenate'};
dat.MW    = [256.42;      254.4;          284.48;     282.46;   280.44;      278.43];
dat.ID    = {'s_3740';    's_3741';       's_3742';   's_3743'; 'm_0128';    'm_0129'};

dat.compId   = {'c','ce','e','er','erm','g','gm','lp','m','mm','n','p','v','vm'};
dat.compName = {'cytoplasm','cell envelope','extracellular','endoplasmic reticulum',...
    'endoplasmic reticulum membrane','Golgi','Golgi membrane','lipid particle','mitochondrion',...
    'mitochondrial membrane','nucleus','peroxisome','vacuole','vacuolar membrane'};

for i=1:length(template.metName)
    chainList = template.chains(i,:);
    chainList = chainList(~cellfun('isempty',chainList));
    for j=1:length(chainList)
        clear rxnsToAdd
        chain = split(chainList(j),',');
        lipid = template.metName(i);
        if startsWith(lipid,'fatty acid')
            idx = find(contains(dat.chain,chain{1}));
            lipid = dat.name(idx);
        else
            for k=1:length(chain)
                lipid = regexprep(lipid,['CHAIN' num2str(k)],chain{k});
            end
        end
        [chainCount,~,n]    = unique(chain);
        [n,~]       = histc(n,unique(n));
        [~,coeff]   = ismember(chainCount,dat.chain);
        
        products    = [template.bbID(i); dat.ID(coeff)];
        totalMW     = sum([template.bbMW(i); dat.MW(coeff).*n]);
        coeff       = [totalMW; dat.MW(coeff).*n]/1000;
        substrate   = model.mets(getIndexes(model,[lipid{1} '[' template.comps{i} ']'], 'metcomps'));
        % Add reactions
        rxnsToAdd.mets          = [substrate; products];
        rxnsToAdd.stoichCoeffs  = [-1; coeff];
        compIdx = find(contains(dat.compId,template.comps{i}));
        rxnsToAdd.rxnNames      = {[lipid{1} ' [' dat.compName{compIdx} '] SLIME rxn']};
        rxnsToAdd.rxnMiriams    = struct('name',{{'sbo'}},'value',{{'SBO:0000395'}});
        rxnsToAdd.lb            = 0;
        rxnsToAdd.rxnConfidenceScores = 1;
        rxnsToAdd.rxns          = generateNewIds(model,'rxns','t_',1);
        if exist('sourcemodel')
            [Lia,Locb]=ismember(rxnsToAdd.rxnNames, sourcemodel.rxnNames);
            if any(Locb)
                rxnsToAdd.rxns(Lia) = sourcemodel.rxns(Locb(Lia));
            end
        end
        model                   = addRxns(model,rxnsToAdd);
    end
end
