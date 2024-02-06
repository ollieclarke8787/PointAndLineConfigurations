load "gridGens.m2"
load "gridExample.m2"

L = gens P2';
M = for j from 0 to 43 list det L^{0}_{j};
-- check that each generator is required 
for i from 0 to 43 do (
    m := M#i;
    N := drop(M, i);
    I := ideal N;
    print not (m%I == 0);
    )
    
