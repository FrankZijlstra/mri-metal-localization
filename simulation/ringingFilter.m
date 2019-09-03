function [ F ] = ringingFilter (s)

% F = 0.75 + 0.25 * repmat(hann(s(1)),[1 s(2) s(3)]) .* repmat(reshape(hann(s(2)),[1 s(2)]),[s(1) 1 s(3)]) .* repmat(reshape(hann(s(3)),[1 1 s(3)]),[s(1) s(2) 1]);

if (length(s) == 3)
    F = repmat(tukeywin(s(1)),[1 s(2) s(3)]) .* repmat(reshape(tukeywin(s(2)),[1 s(2)]),[s(1) 1 s(3)]) .* repmat(reshape(tukeywin(s(3)),[1 1 s(3)]),[s(1) s(2) 1]);
elseif (length(s) == 2)
    F = repmat(tukeywin(s(1)),[1 s(2)]) .* repmat(reshape(tukeywin(s(2)),[1 s(2)]),[s(1) 1]);
end

end
