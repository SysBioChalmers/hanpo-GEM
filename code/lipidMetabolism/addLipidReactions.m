function model = addLipidReactions(template, model, sourcemodel)
% addLipidReactions
%  Uses template reactions in a defined format to add reactions in lipid
%  metabolism, with variable acyl-chain and localization.
%
%  Input:
%    template    structure with defined fields:
%     .rxns      strings with reaction names, with placeholders to be
%                replaced with specified acyl chains and compartments.
%     .eqns      strings with written out reaction equations, with
%                placeholders for acyl-chains and compartments.
%     .comps     strings with metabolite abbreviation(s), separated by
%                comma, as replacement for placeholder.
%     .grRules   strings with grRules, separated by comma, with the
%                same length as the comps field (each localization has
%                its own grRule specified. (optional)
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

dat.chain = {'16:0';          '16:1';             '18:0';         '18:1';       '18:2';          '18:3'};
dat.ate   = {'palmitate';     'palmitoleate';     'stearate';     'oleate';     'linoleate';     'linolenate'};
dat.coa   = {'palmitoyl-CoA'; 'palmitoleoyl-CoA'; 'stearoyl-CoA'; 'oleoyl-CoA'; 'linoleoyl-CoA'; 'linolenoyl-CoA'};

dat.compId   = {'c','ce','e','er','erm','g','gm','lp','m','mm','n','p','v','vm'};
dat.compName = {'cytoplasm','cell envelope','extracellular','endoplasmic reticulum',...
    'endoplasmic reticulum membrane','Golgi','Golgi membrane','lipid particle','mitochondrion',...
    'mitochondrial membrane','nucleus','peroxisome','vacuole','vacuolar membrane'};

for l = 1:numel(template.rxns)
    chains  = template.chains(l,:);
    chains(cellfun('isempty',chains)) = [];
    newEqns = cellstr(repmat(template.eqns{l}, length(chains),1));
    newRxns = cellstr(repmat(template.rxns{l}, length(chains),1));
    for i = 1:numel(newEqns)
        chainComb = strtrim(split(chains{i},','));
        for j = 1:numel(chainComb)
            POS = ['CHAIN' num2str(j)];
            [~, chainLoc]   = ismember(chainComb{j},dat.chain);
            chainAtes       = dat.ate(chainLoc);
            chainCoas       = dat.coa(chainLoc);
            
            newEqns(i) = regexprep(newEqns(i),[POS '_NUMB'],chainComb{j});
            newEqns(i) = regexprep(newEqns(i),[POS '_ATE'],chainAtes);
            newEqns(i) = regexprep(newEqns(i),[POS '_COA'],chainCoas);
            newRxns(i)= regexprep(newRxns(i),[POS '_COA'],chainCoas);
            newRxns(i)= regexprep(newRxns(i),[POS '_ATE'],chainAtes);
            newRxns(i)= regexprep(newRxns(i),POS,chainComb{j});
        end
    end
    if isfield(template,'comps')
        comps = template.comps{l};
        comps = strtrim(split(comps,','));
        [newEqns, newRxns] = fillComps(newEqns, newRxns, comps,dat);
    end
    if isfield(template,'grRules')
        grRules = template.grRules{l};
        grRules = strtrim(split(grRules,','));
        if length(grRules)==1
            newGrRules = cell(repmat(grRules,length(newEqns),1));
        else
            idx = 1:1:length(newEqns);
            idx = reshape(idx,[length(idx)/length(comps),length(comps)]);
            for j = 1:length(comps)
                newGrRules(idx(:,j),1) = grRules(j);
            end
        end
    end
    
    % Query reaction name to see if it already exists in the model. If so,
    % then don't add it.
    [Lia, ~]             = ismember(newRxns, model.rxnNames);
    if isfield(template,'grRules')
        rxnsToAdd.grRules    = newGrRules(~Lia);
    end
    rxnsToAdd.equations  = newEqns(~Lia);
    rxnsToAdd.rxnNames   = newRxns(~Lia);
    rxnsToAdd.rxns       = generateNewIds(model,'rxns','t_',length(rxnsToAdd.equations));
    if exist('sourcemodel')
        [Lia,Locb]=ismember(rxnsToAdd.rxnNames, sourcemodel.rxnNames);
        if any(Locb)
            rxnsToAdd.rxns(Lia) = sourcemodel.rxns(Locb(Lia));
        end
    end
    model                   = addRxns(model,rxnsToAdd,3,'','m_',true);
end
end

%% Duplicate for compartments
function [newEqns, newRxns] = fillComps(templEqns, templNames, comps, dat)

if ~iscell(templEqns)
    templEqns   = {templEqns};
end

[~, compLoc]    = ismember(comps,dat.compId);
compNames       = dat.compName(compLoc);

if ~iscell(comps)
    comps       = {comps};
    newEqns     = cell(length(templEqns),1);
    newEqns(:)  = templEqns;
    newRxns    = cell(length(templNames),1);
    newRxns(:) = templNames;
else
    newEqns     = repmat(templEqns,length(comps),1);
    newRxns    = repmat(templNames,length(comps),1);
end

idx = 1:1:length(newEqns);
idx = reshape(idx,[length(idx)/length(comps),length(comps)]);

for i=1:length(comps)
    newEqns(idx(:,i))  = regexprep(newEqns(idx(:,i)),'\[COMP\]',['\[' comps{i} '\]']);
    newRxns(idx(:,i)) = regexprep(newRxns(idx(:,i)),'COMP',['' compNames{i} '']);
end
end
