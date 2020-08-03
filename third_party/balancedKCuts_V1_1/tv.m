function tot_var = tv(W, f)
%
% (C)2014-15 Syama Sundar Rangapuram, Pramod Kaushik Mudrakarta and Matthias Hein 
% Machine Learning Group, Saarland University, Saarbruecken, Germany
% http://www.ml.uni-saarland.de
%
    
    Wtriu = triu(W);
    [ixt, jxt, wvalt] = find(Wtriu);

    tot_var = wvalt'*abs( f(ixt) - f(jxt) );
    
    if isempty(tot_var)        
        tot_var = 0;
    end
    
end