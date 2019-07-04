%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [model,k] = scaleLipids(model,scaling)
%
%   Modified from SLIMEr, under MIT License:
%   https://github.com/SysBioChalmers/SLIMEr/blob/master/models/scaleAbundancesInModel.m
%
%   Specify in scaling whether 'tails' or 'chains' should be scaled.
%
% 2018-09-25    Eduard Kerkhoven
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [model,k] = scaleLipids(model,scaling)
%Find optimal scaling factor:
k0   = 1;
kOpt = fminsearch(@(k)unusedLipid(k,model,scaling),k0);

%Find optimality range:
krange(1) = fminsearch(@(k) +minScaling(k,model,scaling,kOpt),kOpt);
krange(2) = fminsearch(@(k) -minScaling(k,model,scaling,kOpt),kOpt);
disp(['Optimality range: k = [ ' num2str(krange(1)) ' , ' num2str(krange(2)) ' ]'])
if(krange(1) == krange(2))
    error('Could not find an optimality range!')
end

%Scale with the average of the range:
k     = mean(krange);
model = adjustModel(model,k,true,scaling);
disp(['Scaled ' scaling(1:end-1) ' data in model: k = ' num2str(k)])

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function exchange = unusedLipid(k,model,scaling)

%Adjust stoich coeffs of the corresponding pseudo-rxn:
model = adjustModel(model,k,false,scaling);

%Optimize model:
try
    [sol,~] = solveModel(model);
    if isempty(sol.x)
        sol.x = ones(size(model.rxns));
    end
catch
    sol.x = ones(size(model.rxns));
end

%Objective function: unused tails or backbones
exchange_tails = sol.x(getIndexes(model,'r_4064','rxns'));
exchange_backs = sol.x(getIndexes(model,'r_4062','rxns'));
exchange       = exchange_tails + exchange_backs;

disp(['Scaling abundance data: k = ' num2str(k) ' -> exchange = ' num2str(exchange)])

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function k = minScaling(k,model,scaling,kOpt)

%Adjust stoich coeffs of the corresponding pseudo-rxn:
model = adjustModel(model,k,true,scaling);

%Optimize model:
try
    [sol,~] = solveModel(model);
    ngamIdx = getIndexes(model,'r_4046','rxns');
    disp(['Finding scaling range: k = ' num2str(k) ' -> Maintenance = ' num2str(sol.x(ngamIdx))])
catch
    disp(['Finding scaling range: k = ' num2str(k) ' -> Maintenance = ' num2str(0)])
    k = kOpt;  %any unfeasible simulation returns the original value
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model = adjustModel(model,k,block,scaling)

%Block exchange of tails and backbones:
if block
    backbExIdx  = getIndexes(model,'r_4062','rxns');
    chainExIdx  = getIndexes(model,'r_4064','rxns');
    model = setParam(model,'eq',[backbExIdx,chainExIdx],0);
end

%Switch what to rescale depending on flag:
switch scaling
    case 'backbones'
        scaleRxn = getIndexes(model,'r_4063','rxns');
    case 'tails'
        scaleRxn = getIndexes(model,'r_4065','rxns');
end

%Find positions:
scaleMets = model.S(:,scaleRxn) < 0;

%Change stoich coeffs:
model.S(scaleMets,scaleRxn) = k*model.S(scaleMets,scaleRxn);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [sol,model] = solveModel(model)
model = setParam(model,'obj','r_2111',1);
sol   = solveLP(model);
model = setParam(model,'eq','r_2111',-sol.f);
model = setParam(model,'obj','r_4046',1);
sol = solveLP(model,1);
end