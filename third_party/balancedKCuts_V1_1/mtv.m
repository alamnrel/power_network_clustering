function tot_var = mtv(W, F)
%
% (C)2014-15 Syama Sundar Rangapuram, Pramod Kaushik Mudrakarta and Matthias Hein 
% Machine Learning Group, Saarland University, Saarbruecken, Germany
% http://www.ml.uni-saarland.de
%

    Wtriu = triu(W);
    [ixt, jxt, wvalt] = find(Wtriu);

    k = size(F, 2);
    tot_var = zeros(k, 1);
    for l=1:k
			if ~isempty(wvalt)
				tot_var(l) = wvalt'*abs( F(ixt, l) - F(jxt, l) );
			end
    end
       
end
