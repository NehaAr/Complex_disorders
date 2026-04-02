source("http://zzlab.net/GAPIT/gapit_functions.txt")
myY <- read.table("http://zzlab.net/GAPIT/data/mdp_traits.txt", head = TRUE)
yGD=read.table(file="http://zzlab.net/GAPIT/data/mdp_numeric.txt",head=T)
myGM=read.table(file="http://zzlab.net/GAPIT/data/mdp_SNP_information.txt",head=T)
GAPIT.Power.compare(
GD=myGD,
GM=myGM,
nrep=10,
h2=0.9,
model =c("GLM","MLM","CMLM","MLMM","FarmCPU","BLINK","SUPER"),
NQTN=5)
