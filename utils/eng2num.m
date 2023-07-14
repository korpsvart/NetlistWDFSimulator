function number = eng2num(stringValue)
%ENG2NUM Convert string with engineering-notation to number
%   Engineering notation comes from SPICE3-standard.

% Original function created by Roman Müller-Hainbach (license below).
% Slightly adjusted for our purposes.

% Copyright (c) 2016, Roman Müller-Hainbach
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.


i = length(stringValue);
while i > 0 && isnan(str2double(stringValue(1:i)))
    i = i - 1;
end
firstNumber = stringValue(1:i);
% Extract first number part

stringValue = stringValue(i+1:end);
j = 0;
while j < length(stringValue) && isstrprop(stringValue(j+1), 'alpha')
    j = j + 1;
end
firstLetters = stringValue(1:j);
% Extract first letters part

stringValue = stringValue(j+1:end);
k = 0;
while k < length(stringValue) && ~isnan(str2double(stringValue(1:k+1)))
    k = k + 1;
end
secondNumber = stringValue(1:k);
% Extract second number part

remainder = stringValue(k+1:end);

if isempty(firstNumber)
    number = NaN;
else
    number = str2double(firstNumber);
end
% Evaluate first number part

if ~isempty(secondNumber)
    if contains(firstNumber,'.') || contains(secondNumber,'.')
        number = NaN;
        return;
    end
    if number >= 0
        sign = 1;
    else
        sign = -1;
    end
    number = number + sign*10^-length(secondNumber)*str2double(secondNumber);
end
% Evaluate second number part


factor = 1;
if length(firstLetters) >= 1
    if length(firstLetters) >= 3
        switch lower(firstLetters(1:3))
            case 'meg'
                factor = 10^6;
            case 'mil'
                factor = 25.4 * 10^-6;
        end
        unitInFirstLetters = factor ~= 1 && length(firstLetters) > 3;
    end
    % Check for 3 letter order of magnitude
    if factor == 1
        switch lower(firstLetters(1))
            case 't'
                factor = 10^12;
            case 'g'
                factor = 10^9;
            case 'k'
                factor = 10^3;
            case 'm'
                factor = 10^-3;
            case 'u'
                factor = 10^-6;
            case 'µ'
                factor = 10^-6;
            case 'µ'
                factor = 10^-6;
            case 'n'
                factor = 10^-9;
            case 'p'
                factor = 10^-12;
            case 'f'
                factor = 10^-15;
        end
        % Check for single letter order of magnitude
        unitInFirstLetters = length(firstLetters) > 1;
    end
else
    unitInFirstLetters = false;
end
% Evaluate factor in first letter part

if ~isempty(remainder) && (any(isstrprop(remainder,'digit')) || unitInFirstLetters)
    number = NaN;
    return;
else
    number = number * factor;
end
% Check remainder