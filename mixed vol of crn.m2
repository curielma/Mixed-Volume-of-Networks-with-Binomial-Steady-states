--Compute the mixed volume of a square network
--GOAL: Compute the mixed volume associated to a polynomial system
--    	  given by a CRN under mass-action kinetics
--INPUT: a chemical reaction network G
--OUTPUT: mixed volume of G
--CREDIT: "Mixed volume of small reaction networks" by
--    	      Nida Obatake, Anne Shiu, and Dilruba Sofia
loadPackage "ReactionNetworks"
loadPackage "LatticePolytopes"

G = reactionNetwork{"A+B --> B+C","B+C --> C+D","C+D --> A+D","A+D --> A+B"}

createRing(G,QQ)

f1 = subRandomReactionRates G
f2 = subRandomInitVals G

R = QQ[G.ConcentrationRates]
f1 = apply(f1, p -> sub(p,R))
f2 = apply(f2, p -> sub(p,R))

W = transpose gens ker transpose stoichiometricMatrix G

(Q,L,U) = LUdecomposition( sub(W,RR) )

f = f1;
for i from 0 to numRows W - 1 do(f = replace( position(flatten entries U^{i}, j -> j != 0), f2_i, f))

mixedVolume apply(f, p -> newtonPolytope p)


****************************************************************************************************
--Network Examples--

--Edelstein
G = reactionNetwork{"A <--> 2A","A+B <--> C","C <--> B"}

--Modified Edelstein
G = reactionNetwork{"A+C <--> 2A+C","A+B <--> C","C <--> B"}
