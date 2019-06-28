equationStrings = {};
%check exchange reactions
    for eqn = 1:length(sol.x(getIndexes(modelHpoGF,sol.x,'rxns')))
        equationStrings{eqn,1}=constructEquations(superModel.subModels{organism},metIdx(upMet,organism));
    end
metEquations = equationStrings;