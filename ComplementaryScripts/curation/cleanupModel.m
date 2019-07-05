% Fix redundancies in complex grRules
for n = 1:length(model.grRules)
    if any(model.grRules{n})
        noAnd = strfind(model.grRules(n),'and');
        noAnd = any(vertcat(noAnd{:})); % Give 0 if no 'and' is present.
        if noAnd == 0
            geneList = transpose(cell(unique(regexp(model.grRules{n},'[)(]*|( and )*|( or )*','split'))));
            geneList = regexprep(geneList,'[(*)*]','');
            if length(geneList) == 1
                newgrRule = geneList;
            else
                newgrRule = geneList{1};
                for k = 2:length(geneList)
                    newgrRule = [newgrRule ' or ' geneList{k}];
                end
            end
            model.grRules(n) = cellstr(newgrRule);
        end
    end
end


%% Remove 'sce' from subsystems
model.subSystems = cellfun(@(x) regexprep(x,'sce[0-9]+ +',''),model.subSystems, 'UniformOutput', 0);

%% Remove unused metabolites
model = removeMets(model,all(model.S == 0,2),false,true,true,true);
