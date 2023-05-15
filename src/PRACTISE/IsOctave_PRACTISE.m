% Haerer, Bernhardt and Schulz (2016)
% "PRACTISE - Photo Rectification And ClassificaTIon SoftwarE (V.2.1)"
%
%   written by
%   Stefan Haerer (LMU Munich)
%   11/2015
%   contact: stefan.haerer@boku.ac.at
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [RetVal] = IsOctave_PRACTISE()
%   Name:       IsOctave_PRACTISE
%   Purpose:    Check if PRACTISE is run in Octave
%
%   Output:     RetVal (Return true if the environment is Octave)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
persistent CacheVal;  % speeds up repeated calls
%
if isempty (CacheVal)
    CacheVal = (exist ('OCTAVE_VERSION', 'builtin') > 0);
end
%
RetVal = CacheVal;
