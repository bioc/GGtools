library(GGtools)
library(ff)

tenOn2021 =
c( "GI_4557248-S",  "GI_15451784-S", "GI_4557678-S", "GI_9951914-S", "GI_21327679-S",
 "GI_7669476-I", "GI_7669478-A", "GI_4557290-I", "GI_4557294-A", "GI_41406053-S")

c22 = getSS("GGdata", "22")
c22 = c22[probeId(tenOn2021),]

e1 = eqtlTests(c22, ~1, targdir = "ooo")
m1 = meqtlTests(list(c22, c22), list(~1, ~1), targdir = "ooobb")

ae = as.ram(e1@fffile)
am = as.ram(m1@fffile)
all.equal(as.numeric(ae)*2 ,as.numeric(am) )

system("rm -rf ooo")
system("rm -rf ooobb")
