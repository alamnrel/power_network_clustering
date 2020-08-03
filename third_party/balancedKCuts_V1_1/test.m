try
    load iris
    clustering = balancedKCut(W, k);
    cnstrClustering = transductiveBalancedKCut(W, k, Labels);
    display('test script ran successfully');
catch
    display('something went wrong... check ReadMe.txt for more information...')
end