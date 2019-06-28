
%Delete model.grRules (redundant and possibly conflicting with model.rules):
if isfield(model,'grRules')
    model = rmfield(model,'grRules');
end


%Check if model is a valid SBML structure:
writeCbModel(model,'sbml','tempModel.xml');
[~,errors] = TranslateSBML('tempModel.xml');
if ~isempty(errors)
    delete('tempModel.xml');
    error('Model should be a valid SBML structure. Please fix all errors before saving.')
end
