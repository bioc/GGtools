
t20 = GGtools:::getCisMap()

NL = GGtools:::namelist(t20)

NL1 = NL[[1]]

t20sl = t20@snplocs[NL1]

probe1 = names(NL)[1]

TARG = t20@generanges[probe1]-50000

d = distance(TARG, t20sl)

all(d <= 50000)
