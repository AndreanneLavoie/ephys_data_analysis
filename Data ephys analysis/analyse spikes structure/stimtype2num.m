function output = stimtype2num(stim_type)
%Convert stim type to number for filtspike function
%to pull out tuning type

%spatial frquency: 'sf' = 1
%temporal frequency: 'tf' = 2
%orientation/derection selectivigty 'or' = 3

switch stim_type
    case {'sf', 'SF'}
        output = 1;
    case {'tf', 'TF'}
        output = 2;
    case {'or', 'OR'}
        output = 3;
end

