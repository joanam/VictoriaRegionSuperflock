
# Libraries required
require("readxl")
require("ape")
require("phylobase")
require("tidyr")
library("ggplot2")
library("ggtree")
library("tidytree")
library("ggtreeExtra")
library("gplots")
library("ggnewscale")
library("phytools")
library("ggpattern")
library("Rcpp")
library("dplyr")

# Read in sample information
sampleInfo<-as.data.frame(read_excel("D:/Dropbox/CichlidGenomesProject/planning/Cichlids_Sequenced_New.xlsx"))
groups<-read.table("D:/Dropbox/victoriaGenomes/phylogenetics/lakes.tree.groups.ecology.ind",header=F)
names(groups)<-c("sample","U","group")

# Read in the tree
tree<-read.tree("D:/Dropbox/victoriaGenomes/phylogenetics/allGenomes.1kbThinned.withSRA.chr1-22.alluaOut.noSquami.phylip.varsites.phy.treefile")

# Root the tree on Astatoreochromis alluaudi
outgroups<-c("103218","105845","106754","131201","10685","10719","103881","131344")
tree<-root(tree,outgroup = outgroups,edgelabel = TRUE)

# Ladderize the tree
tree<-ladderize(tree)

# Reorder the sample information to match the order in the tree labels
treetips<-as.character(tree$tip.label)
treetipsorder=cbind(treetips,"treenames"=treetips,"indexx"=c(100:(length(treetips)+99)))
mergedSampleInfo<-merge(groups,treetipsorder,by.x="sample",by.y="treetips",all.y=T)
sortedSampleInfo<-mergedSampleInfo[order(mergedSampleInfo$indexx),]
sortedSampleInfo<-droplevels(sortedSampleInfo)
sortedSampleInfo<-unique(x=sortedSampleInfo)
sortedSampleInfo<-cbind(sortedSampleInfo,"group_color"=sortedSampleInfo$group)
sortedSampleInfo$group_color<-as.factor(sortedSampleInfo$group_color)
colors<-colorRampPalette(c("black","yellow","red","blue","orange"))
levels(sortedSampleInfo$group_color)<-colors(length(levels(sortedSampleInfo$group_color)))
sortedSampleInfo$group_color<-as.character(sortedSampleInfo$group_color)
sortedSampleInfo$group_pch<-as.factor(sortedSampleInfo$group)
levels(sortedSampleInfo$group_pch)<-rep(0:25,times=100)[1:length(levels(sortedSampleInfo$group_pch))]
sortedSampleInfo$group_pch<-as.integer(sortedSampleInfo$group_pch)
mergedSampleInfo<-merge(sampleInfo,sortedSampleInfo,by="sample",all.y=T)
sortedSampleInfo<-mergedSampleInfo[order(mergedSampleInfo$indexx),]
sortedSampleInfo<-droplevels(sortedSampleInfo)
sortedSampleInfo<-unique(x=sortedSampleInfo)

# Add species symbol
sortedSampleInfo$SpeciesPch<-sortedSampleInfo$Species_Name_CleanDM
sortedSampleInfo$SpeciesPch[(sortedSampleInfo$Group!="LVRS Victoria"&sortedSampleInfo$Group!="LVRS Kyoga" )| is.na(sortedSampleInfo$Major_phylogroup)]<-NA
sortedSampleInfo$SpeciesPch<-as.factor(sortedSampleInfo$SpeciesPch)
levels(sortedSampleInfo$SpeciesPch)<-rep(1:20,length.out=length(levels(sortedSampleInfo$SpeciesPch)))
sortedSampleInfo$SpeciesPch<-as.integer(sortedSampleInfo$SpeciesPch)

# remove ecology information for outgroups
sortedSampleInfo$simplifiedEcology[sortedSampleInfo$Lake=="Malawi"|
                            sortedSampleInfo$Species=="calliptera"|
                            sortedSampleInfo$Species=="burtoni" |
                            sortedSampleInfo$Species=="alluaudi"]<-NA
sortedSampleInfo$simplifiedEcology[sortedSampleInfo$Species=="sp. Red cheek"|
                           sortedSampleInfo$Species=="bloyeti"|
                           sortedSampleInfo$Species=="sparsidens" |
                           sortedSampleInfo$Species=="sp. Wami"]<-NA
sortedSampleInfo$simplifiedEcology[sortedSampleInfo$Species=="stappersi"|
                                     sortedSampleInfo$Species=="sp. Yaekama"|
                                     sortedSampleInfo$Species=="gracilior" |
                                     sortedSampleInfo$Species=="pharyngalis"]<-NA


# Add lake colours
sortedSampleInfo<-cbind(sortedSampleInfo,"lake_color"=sortedSampleInfo$Lake)
sortedSampleInfo$lake_color<-as.factor(sortedSampleInfo$lake_color)
colors<-colorRampPalette(c("black","red","grey","blue","orange"))
levels(sortedSampleInfo$lake_color)[levels(sortedSampleInfo$lake_color)=="Albert"]<-"cornflowerblue"
levels(sortedSampleInfo$lake_color)[levels(sortedSampleInfo$lake_color)=="Edward"]<-"blue"
levels(sortedSampleInfo$lake_color)[levels(sortedSampleInfo$lake_color)=="Egypt"]<-"#40A3DBB3"
levels(sortedSampleInfo$lake_color)[levels(sortedSampleInfo$lake_color)=="Sudan"]<-"#40A3DBB3"
levels(sortedSampleInfo$lake_color)[levels(sortedSampleInfo$lake_color)=="Kivu"]<-"darkblue"
levels(sortedSampleInfo$lake_color)[levels(sortedSampleInfo$lake_color)=="Kivu outflow"]<-"darkblue"
levels(sortedSampleInfo$lake_color)[levels(sortedSampleInfo$lake_color)=="Kagera"]<-"red"
levels(sortedSampleInfo$lake_color)[levels(sortedSampleInfo$lake_color)=="Saka"]<-"deepskyblue"
levels(sortedSampleInfo$lake_color)[levels(sortedSampleInfo$lake_color)=="Victoria"]<-"#FF9F00"
levels(sortedSampleInfo$lake_color)[levels(sortedSampleInfo$lake_color)=="Kanyaboli"]<-"#F4C430"
levels(sortedSampleInfo$lake_color)[levels(sortedSampleInfo$lake_color)=="Kyoga"]<-"#FF6304"
levels(sortedSampleInfo$lake_color)[levels(sortedSampleInfo$lake_color)=="Congo"]<-"darkred"
levels(sortedSampleInfo$lake_color)[levels(sortedSampleInfo$lake_color)=="Mweru"]<-"black"
levels(sortedSampleInfo$lake_color)[levels(sortedSampleInfo$lake_color)=="Eastern"]<-"black"
levels(sortedSampleInfo$lake_color)[levels(sortedSampleInfo$lake_color)=="Cunene River"]<-"black"
levels(sortedSampleInfo$lake_color)[levels(sortedSampleInfo$lake_color)=="ref"]<-"black"
levels(sortedSampleInfo$lake_color)[levels(sortedSampleInfo$lake_color)=="Manyara"]<-"black"
levels(sortedSampleInfo$lake_color)[levels(sortedSampleInfo$lake_color)=="Burigi"]<-"#F4C430"
levels(sortedSampleInfo$lake_color)[levels(sortedSampleInfo$lake_color)=="Ikimba"]<-"#F4C430"
levels(sortedSampleInfo$lake_color)[levels(sortedSampleInfo$lake_color)=="Nabugabo"]<-"#F4C430"
levels(sortedSampleInfo$lake_color)[levels(sortedSampleInfo$lake_color)=="Turkana"]<-"#40A3DBB3"
levels(sortedSampleInfo$lake_color)[levels(sortedSampleInfo$lake_color)=="Boukou"]<-"#40A3DBB3"
levels(sortedSampleInfo$lake_color)[levels(sortedSampleInfo$lake_color)=="Mpanga River"]<-"deepskyblue"
levels(sortedSampleInfo$lake_color)[levels(sortedSampleInfo$lake_color)=="Victoria-allua"]<-"black"
levels(sortedSampleInfo$lake_color)[levels(sortedSampleInfo$lake_color)=="Edward-Nile"]<-"darkblue"
levels(sortedSampleInfo$lake_color)[levels(sortedSampleInfo$lake_color)=="Kivu-Nile"]<-"darkblue"
levels(sortedSampleInfo$lake_color)[levels(sortedSampleInfo$lake_color)=="Kivu-Nile?"]<-"darkblue"
levels(sortedSampleInfo$lake_color)[levels(sortedSampleInfo$lake_color)=="Victoria-Pseudo"]<-"black"
levels(sortedSampleInfo$lake_color)[levels(sortedSampleInfo$lake_color)=="Kagera-burt"]<-"black"
levels(sortedSampleInfo$lake_color)[levels(sortedSampleInfo$lake_color)=="Kagera-allua"]<-"black"
levels(sortedSampleInfo$lake_color)[levels(sortedSampleInfo$lake_color)=="Boukou-Nile"]<-"darkblue"
levels(sortedSampleInfo$lake_color)[levels(sortedSampleInfo$lake_color)=="Rukwa"]<-"black"
levels(sortedSampleInfo$lake_color)[levels(sortedSampleInfo$lake_color)=="Burungi"]<-"black"
levels(sortedSampleInfo$lake_color)[levels(sortedSampleInfo$lake_color)=="Malawi"]<-"black"
levels(sortedSampleInfo$lake_color)[levels(sortedSampleInfo$lake_color)=="Kumba"]<-"black"
levels(sortedSampleInfo$lake_color)[levels(sortedSampleInfo$lake_color)==""]<-"black"
sortedSampleInfo$lake_color<-as.character(sortedSampleInfo$lake_color)
sortedSampleInfo$lake_color[sortedSampleInfo$Group=="Astatotilapia (Nile)"]<-"black"
sortedSampleInfo$lake_color[sortedSampleInfo$Species=="burtoni"]<-"black"
sortedSampleInfo$lake_color[sortedSampleInfo$Location=="Sweya"]<-"#F4C430"
sortedSampleInfo$lake_color[sortedSampleInfo$Species=="alluaudi"]<-"black"
sortedSampleInfo$lake_color[sortedSampleInfo$Species=="sp. Red cheek"]<-"black"
sortedSampleInfo[sortedSampleInfo$simplifiedEcology=="?"&!is.na(sortedSampleInfo$simplifiedEcology),"simplifiedEcology"]<-NA

# Add genus colour
sortedSampleInfo$GenusCol<-sortedSampleInfo$Genus
sortedSampleInfo$GenusCol<-as.factor(sortedSampleInfo$GenusCol)
levels(sortedSampleInfo$GenusCol)[grep(levels(sortedSampleInfo$GenusCol),pattern="Pundamilia")]<-"Pundamilia"
levels(sortedSampleInfo$GenusCol)[grep(levels(sortedSampleInfo$GenusCol),pattern="Yssichromis")]<-"Yssichromis"
levels(sortedSampleInfo$GenusCol)[grep(levels(sortedSampleInfo$GenusCol),pattern="Enterochromis")]<-"Enterochromis"
levels(sortedSampleInfo$GenusCol)[grep(levels(sortedSampleInfo$GenusCol),pattern="Paralabidochromis")]<-"Paralabidochromis"
sortedSampleInfo$GenusSimplified<-sortedSampleInfo$GenusCol

sortedSampleInfo$GenusSimplified<-as.character(sortedSampleInfo$GenusSimplified)
sortedSampleInfo$GenusSimplified[grepl("Macropleurodus.*Hoplotilapia)",sortedSampleInfo$GenusSimplified)]<-"Hoplotilapia"
sortedSampleInfo$GenusSimplified[grepl("Macropleurodus.*Ptyochromis)",sortedSampleInfo$GenusSimplified)]<-"Ptyochromis"
sortedSampleInfo$GenusSimplified[grepl("Macropleurodus.*Macropleurodus)",sortedSampleInfo$GenusSimplified)]<-"Macropleurodus"
sortedSampleInfo$GenusSimplified[sortedSampleInfo$GenusSimplified=="\"Astatotilapia\""]<-"Astatotilapia"
sortedSampleInfo$GenusSimplified[sortedSampleInfo$GenusSimplified=="Harpago/Prognathochromis"]<-"Harpagochromis"
sortedSampleInfo$GenusSimplified[grepl("\"Thoracochromis\"?",sortedSampleInfo$GenusSimplified)]<-"Thoracochromis"
sortedSampleInfo$GenusSimplified<-as.factor(sortedSampleInfo$GenusSimplified)


sortedSampleInfo$genusCol<-sortedSampleInfo$GenusSimplified
sortedSampleInfo$genusCol<-as.character(sortedSampleInfo$genusCol)
sortedSampleInfo$genusCol[sortedSampleInfo$Lake=="Malawi"]<-NA
sortedSampleInfo$genusCol[sortedSampleInfo$Species=="calliptera"]<-NA
sortedSampleInfo$genusCol[sortedSampleInfo$Species=="burtoni"]<-NA
sortedSampleInfo$genusCol[sortedSampleInfo$genusCol=="Incertae sedis"]<-NA
sortedSampleInfo$genusCol[sortedSampleInfo$genusCol=="\"Neochromis\""]<-"Neochromis"
sortedSampleInfo$genusCol[sortedSampleInfo$Species=="alluaudi"]<-NA
sortedSampleInfo$genusCol[sortedSampleInfo$Species=="stappersi"]<-NA
sortedSampleInfo$genusCol[sortedSampleInfo$Species=="alluaudi"]<-NA

sortedSampleInfo$Lake[is.na(sortedSampleInfo$Lake)]<-""


# Check if the species are monophyletic
speciesInfo<-cbind("sample"=tree$tip.label,
                   "Genus"=as.character(sortedSampleInfo$GenusSimplified),
                   "Species"=sortedSampleInfo$Species_Name_CleanDM,
                   "SpeciesFull"=sortedSampleInfo$Species_Name_CleanDM,
                   "Group"=sortedSampleInfo$Group,
                   "Lake"=sortedSampleInfo$Lake,
                   "col"=sortedSampleInfo$lake_color,
                   "GenusCol"=sortedSampleInfo$genusCol,
                   "EcoCol"=sortedSampleInfo$simplifiedEcology)
speciesInfo<-as.data.frame(speciesInfo)

# replace NA groups by empty entry:
speciesInfo[is.na(speciesInfo$Group),"Group"]<-""

# set species of lake radiations to lake name
speciesInfo$Species[!speciesInfo$Group%in%c("LVRS Victoria","LVRS Kyoga")]<-
  speciesInfo$Lake[!speciesInfo$Group%in%c("LVRS Victoria","LVRS Kyoga")]

speciesInfo$Species[speciesInfo$Lake%in%c("Ikimba","Burigi")]<-
  speciesInfo$Lake[speciesInfo$Lake%in%c("Ikimba","Burigi")]
speciesInfo<-speciesInfo[speciesInfo$Species!="egg sneaker",]
speciesInfo[is.na(speciesInfo$Species),"Species"]<-
  sortedSampleInfo[is.na(speciesInfo$Species),"group"]
speciesInfo$Species[speciesInfo$sample=="ERR3634113"]<-"Rufiji"
speciesInfo$Species[speciesInfo$sample=="ERR715538"]<-"Malawi"
speciesInfo$Species[speciesInfo$Genus=="Astatoreochromis"]<-"Astatoreochromis"
speciesInfo$Species[speciesInfo$Species=="Lipochromis sp. \"velvet black cryptodon\""]<-"Lipochromis cryptodon"
speciesInfo$Species[speciesInfo$Species=="Mpanga River"]<-"Saka"
speciesInfo[speciesInfo$sample=="13653","Species"]<-"Pundamilia macrocephala Makobe"
speciesInfo[speciesInfo$Species=="Pundamilia macrocephala","Species"]<-"Pundamilia macrocephala Python"
speciesInfo[speciesInfo$sample=="131282","Species"]<-"Astatotilapia burtoni"
speciesInfo[speciesInfo$sample=="Aburtoni","Species"]<-"Astatotilapia burtoni"
speciesInfo[speciesInfo$sample=="Aburtoni","Species"]<-"Astatotilapia burtoni"
speciesInfo[sortedSampleInfo$Group=="Astatotilapia (Nile)"&!is.na(sortedSampleInfo$Group),"Species"]<-"Upper Nile A gracilior"
speciesInfo[sortedSampleInfo$Species=="pharyngalis","Species"]<-"Upper Nile T pharyngalis"
speciesInfo[speciesInfo$sample=="ERR715501","Species"]<-"Ruaha"
speciesInfo[speciesInfo$Species=="Kumba","Species"]<-"Eastern"
speciesInfo[speciesInfo$Species=="Burungi","Species"]<-"Eastern"
speciesInfo[speciesInfo$Species=="Manyara","Species"]<-"Eastern"
speciesInfo[speciesInfo$sample=="81309","Species"]<-"Congo Yaekama"
speciesInfo[speciesInfo$Species=="Congo","Species"]<-"Congo A stappersi"
speciesInfo[speciesInfo$sample=="80054","Species"]<-speciesInfo[speciesInfo$sample=="104063","Species"]
speciesInfo[speciesInfo$sample=="130859","Species"]<-"Boukou red"
speciesInfo[speciesInfo$sample=="130858","Species"]<-"Boukou blue"
speciesInfo[speciesInfo$sample=="103767","Species"]<-"Enterochromis coprologus II"
speciesInfo$Species[grepl(speciesInfo$Species,pattern="Mpan")]<-"Saka"


# Select some Kivu species:
speciesInfo[speciesInfo$Group=="LVRS Kivu"&speciesInfo$Genus=="\"Neochromis\"","Species"]<-"\"Neo\" olivaceus Kivu"
speciesInfo[speciesInfo$Group=="LVRS Kivu"&speciesInfo$Genus=="\"Prognathochromis\"","Species"]<-"Pro vittatus Kivu"
speciesInfo[speciesInfo$Group=="LVRS Kivu"&speciesInfo$Genus=="\"Lipochromis\"","Species"]<-"Lip occultidens Kivu"
speciesInfo[speciesInfo$Group=="LVRS Kivu"&grepl(speciesInfo$SpeciesFull,pattern="Pundamilia-like sp1"),"Species"]<-"inc littoral black sp. 1 Kivu"
speciesInfo[speciesInfo$Group=="LVRS Kivu"&grepl(speciesInfo$SpeciesFull,pattern="graueri"),"Species"]<-"Psa graueri Kivu"
speciesInfo[speciesInfo$Group=="LVRS Kivu"&grepl(speciesInfo$SpeciesFull,pattern="kami"),"Species"]<-"Yss kamiranzovu Kivu"
speciesInfo[speciesInfo$Group=="LVRS Kivu"&grepl(speciesInfo$SpeciesFull,pattern="pauc"),"Species"]<-"Par paucidens Kivu"

# Select some Edward species
speciesInfo[speciesInfo$Group=="LVRS Edward"&grepl(speciesInfo$SpeciesFull,pattern="Psamm"),"Species"]<-"Psa schubotzi complex Edward"
speciesInfo[speciesInfo$Group=="LVRS Edward"&grepl(speciesInfo$SpeciesFull,pattern="oregosoma"),"Species"]<-"Ast oregosoma Edward"
speciesInfo[speciesInfo$Group=="LVRS Edward"&grepl(speciesInfo$SpeciesFull,pattern="Lipo"),"Species"]<-"Lip taurinus Edward"
speciesInfo[speciesInfo$Group=="LVRS Edward"&grepl(speciesInfo$SpeciesFull,pattern="limax"),"Species"]<-"Hap limax Edward"

# Select some Kagera species
speciesInfo[speciesInfo$Group=="LVRS Kagera"&grepl(speciesInfo$SpeciesFull,pattern="Progn"),"Species"]<-"Prognathochromis Kagera"
speciesInfo[speciesInfo$Group=="LVRS Kagera"&grepl(speciesInfo$SpeciesFull,pattern="Gau"),"Species"]<-"Gau sp. \"red top\" Kagera"
speciesInfo[speciesInfo$Group=="LVRS Kagera"&grepl(speciesInfo$SpeciesFull,pattern="Ente"),"Species"]<-"Ent sp. \"Ngoma\" Kagera"
speciesInfo[speciesInfo$Group=="LVRS Kagera"&grepl(speciesInfo$SpeciesFull,pattern="Psa"),"Species"]<-"Psa sp. \"red\" Kagera"
speciesInfo[speciesInfo$Group=="LVRS Kagera"&grepl(speciesInfo$SpeciesFull,pattern="Merule"),"Species"]<-"Ast spp. Merule"

# Albert species
speciesInfo[speciesInfo$Group=="LVRS Albert"&grepl(speciesInfo$SpeciesFull,pattern="aeneocolor-like 1"),"Species"]<-"Ast sp. \"aeneocolor-like 1\" Albert"
speciesInfo[speciesInfo$Group=="LVRS Albert"&grepl(speciesInfo$SpeciesFull,pattern="Neochromis "),"Species"]<-"\"Neo\" sp. Albert"
speciesInfo[speciesInfo$Group=="LVRS Albert"&grepl(speciesInfo$SpeciesFull,pattern="avium"),"Species"]<-"Tho avium Albert"
speciesInfo[speciesInfo$Group=="LVRS Albert"&grepl(speciesInfo$SpeciesFull,pattern="wingatii"),"Species"]<-"Tho wingatii Albert"
speciesInfo[speciesInfo$Group=="LVRS Albert"&grepl(speciesInfo$SpeciesFull,pattern="mahagiensis"),"Species"]<-"Tho mahagiensis Albert"
speciesInfo[speciesInfo$Group=="LVRS Albert"&grepl(speciesInfo$SpeciesFull,pattern="ocolor-like 2"),"Species"]<-"Ast sp. \"aeneocolor-like 2\" Albert"


# Simplify labels
tree$tip.label<-paste(sortedSampleInfo$sample,sortedSampleInfo$GenusSimplified,sortedSampleInfo$Species)

sortedSampleInfo<-cbind(tiplab=tree$tip.label,sortedSampleInfo)
sortedSampleInfo$SpShape<-as.character(sortedSampleInfo$SpeciesPch)

sortedSampleInfo$spLabel<-paste(sortedSampleInfo$Species_Name_CleanDM)

# Abbreviate the genera
sortedSampleInfo$spLabel<-gsub('Pundamilia', 'Pun', sortedSampleInfo$spLabel)
sortedSampleInfo$spLabel<-gsub('Neochromis', 'Neo', sortedSampleInfo$spLabel)
sortedSampleInfo$spLabel<-gsub('Yssichromis', 'Yss', sortedSampleInfo$spLabel)
sortedSampleInfo$spLabel<-gsub('Prognathochromis', 'Pro', sortedSampleInfo$spLabel)
sortedSampleInfo$spLabel<-gsub('Harpagochromis', 'Har', sortedSampleInfo$spLabel)
sortedSampleInfo$spLabel<-gsub('Tridontochromis', 'Tri', sortedSampleInfo$spLabel)
sortedSampleInfo$spLabel<-gsub('Lithochromis', 'Lit', sortedSampleInfo$spLabel)
sortedSampleInfo$spLabel<-gsub('Mbipia', 'Mbi', sortedSampleInfo$spLabel)
sortedSampleInfo$spLabel<-gsub('Paralabidochromis', 'Par', sortedSampleInfo$spLabel)
sortedSampleInfo$spLabel<-gsub('Incertae sedis', 'inc', sortedSampleInfo$spLabel)
sortedSampleInfo$spLabel<-gsub('Mac super−lineage (Pty)', 'Pty', sortedSampleInfo$spLabel)
sortedSampleInfo$spLabel<-gsub('Mac super−lineage (Hop)', 'Hop', sortedSampleInfo$spLabel)
sortedSampleInfo$spLabel<-gsub('Mac super−lineage (Mac)', 'Mac', sortedSampleInfo$spLabel)
sortedSampleInfo$spLabel<-gsub('Ptyochromis', 'Pty', sortedSampleInfo$spLabel)
sortedSampleInfo$spLabel<-gsub('Platytaeniodus', 'Pla', sortedSampleInfo$spLabel)
sortedSampleInfo$spLabel<-gsub('Astatotilapia', 'Ast', sortedSampleInfo$spLabel)
sortedSampleInfo$spLabel<-gsub('Enterochromis', 'Ent', sortedSampleInfo$spLabel)
sortedSampleInfo$spLabel<-gsub('Hoplotilapia', 'Hop', sortedSampleInfo$spLabel)
sortedSampleInfo$spLabel<-gsub('Gaurochromis', 'Gau', sortedSampleInfo$spLabel)
sortedSampleInfo$spLabel<-gsub('Labrochromis', 'Lab', sortedSampleInfo$spLabel)
sortedSampleInfo$spLabel<-gsub('Lipochromis', 'Lip', sortedSampleInfo$spLabel)
sortedSampleInfo$spLabel<-gsub('Macropleurodus', 'Mac', sortedSampleInfo$spLabel)
sortedSampleInfo$spLabel<-gsub('Psammochromis', 'Psa', sortedSampleInfo$spLabel)
sortedSampleInfo$spLabel<-gsub('Double-stripe group', 'dou', sortedSampleInfo$spLabel)
sortedSampleInfo$spLabel<-gsub('Pyxichromis', 'Pyx', sortedSampleInfo$spLabel)
sortedSampleInfo$spLabel<-gsub('Haplochromis', 'Hap', sortedSampleInfo$spLabel)
sortedSampleInfo$spLabel<-gsub('Thoracochromis', 'Tho', sortedSampleInfo$spLabel)
sortedSampleInfo$spLabel[!grepl(sortedSampleInfo$Lake,pattern="Victoria|Kyoga")]<-
  paste(sortedSampleInfo$spLabel,sortedSampleInfo$Lake)[!grepl(sortedSampleInfo$Lake,pattern="Victoria|Kyoga")]
sortedSampleInfo$spLabel[sortedSampleInfo$Lake=="Malawi"]<-NA
sortedSampleInfo$spLabel[sortedSampleInfo$Species=="burtoni"]<-NA
sortedSampleInfo$spLabel[sortedSampleInfo$Species=="calliptera"]<-NA
sortedSampleInfo$spLabel[sortedSampleInfo$Species=="alluaudi"]<-NA
sortedSampleInfo$spLabel[sortedSampleInfo$Lake=="Turkana"&!is.na(sortedSampleInfo$Lake)]<-"Lake Turkana"
sortedSampleInfo$spLabel[sortedSampleInfo$Lake=="Egypt"&!is.na(sortedSampleInfo$Lake)]<-"Egypt"
sortedSampleInfo$spLabel[sortedSampleInfo$Lake=="Sudan"&!is.na(sortedSampleInfo$Lake)]<-"Sudan"
sortedSampleInfo$spLabel[sortedSampleInfo$Lake=="Mpanga River"&!is.na(sortedSampleInfo$Lake)]<-"Lake Saka / Mpanga River"
sortedSampleInfo$spLabel[sortedSampleInfo$Lake=="Saka"&!is.na(sortedSampleInfo$Lake)]<-"Lake Saka / Mpanga River"
sortedSampleInfo$spLabel[sortedSampleInfo$Lake=="Albert"&!is.na(sortedSampleInfo$Lake)]<-"Lake Albert"
sortedSampleInfo$spLabel[sortedSampleInfo$Group=="LVRS Edward"&!is.na(sortedSampleInfo$Group)]<-"Lake Edward"
sortedSampleInfo$spLabel[sortedSampleInfo$Group=="LVRS Kivu"&!is.na(sortedSampleInfo$Group)]<-"Lake Kivu"
sortedSampleInfo$spLabel[sortedSampleInfo$Group=="LVRS Kagera"&!is.na(sortedSampleInfo$Group)]<-"Kagera Lakes"
sortedSampleInfo$spLabel[sortedSampleInfo$Lake=="Burigi"&!is.na(sortedSampleInfo$Lake)]<-"Ast nubila Burigi"
sortedSampleInfo$spLabel[sortedSampleInfo$Lake=="Ikimba"&!is.na(sortedSampleInfo$Lake)]<-"Ast nubila Ikimba"
sortedSampleInfo$spLabel[sortedSampleInfo$Lake=="Rukwa"&!is.na(sortedSampleInfo$Lake)]<-"Lake Rukwa"

# Add label line colour
sortedSampleInfo$labelCol<-ifelse(is.na(sortedSampleInfo$spLabel),"white","grey")
sortedSampleInfo$labelCol[sortedSampleInfo$spLabel=="white"]<-"white"

# Add species colour
sortedSampleInfo$SpeciesCol<-sortedSampleInfo$spLabel
sortedSampleInfo$SpeciesCol<-as.factor(sortedSampleInfo$SpeciesCol)
levels(sortedSampleInfo$SpeciesCol)[levels(sortedSampleInfo$SpeciesCol)%in%names(table(sortedSampleInfo$SpeciesCol)[table(sortedSampleInfo$SpeciesCol)==1])]<-NA
levels(sortedSampleInfo$SpeciesCol)[!grepl(levels(sortedSampleInfo$SpeciesCol),pattern="grey")]<-colors(length(levels(sortedSampleInfo$SpeciesCol)[!grepl(levels(sortedSampleInfo$SpeciesCol),pattern="grey")]))
sortedSampleInfo$SpeciesCol<-as.character(sortedSampleInfo$SpeciesCol)
sortedSampleInfo$SpeciesCol[sortedSampleInfo$SpeciesCol=="grey"]<-NA
# sortedSampleInfo$SpeciesCol[sortedSampleInfo$Lake=="Albert"]<-"cornflowerblue"
# sortedSampleInfo$SpeciesCol[sortedSampleInfo$Lake=="Edward"]<-"blue"
# sortedSampleInfo$SpeciesCol[sortedSampleInfo$Lake=="Kivu"]<-"darkblue"
# sortedSampleInfo$SpeciesCol[sortedSampleInfo$Lake=="Saka"]<-"deepskyblue"
# sortedSampleInfo$SpeciesCol[sortedSampleInfo$Lake=="Mpanga River"]<-"deepskyblue"
# sortedSampleInfo$SpeciesCol[sortedSampleInfo$Lake=="Rukwa"]<-"black"
sortedSampleInfo$SpeciesCol[!(sortedSampleInfo$Lake=="Victoria"|sortedSampleInfo$Lake=="Kyoga")]<-NA

# Keep only one label per species
# indices for labels
middle<-function(x){x[as.integer(length(x)/2+0.5)]}
a<-numeric()
d<-0
specList<-unique(sortedSampleInfo$spLabel[!is.na(sortedSampleInfo$spLabel)])
for(species in specList){
  x<-which(sortedSampleInfo$spLabel==species)
  d=d+1
  a<-c(a,middle(x))
}
indices<-1:length(tree$tip.label)
toRemove<-indices[-c(a)]

sortedSampleInfo$spLabel[toRemove]<-NA


# Remove species colours for outgroups
sortedSampleInfo$SpeciesCol[sortedSampleInfo$Lake=="Malawi"]<-NA
sortedSampleInfo$SpeciesCol[sortedSampleInfo$Species=="burtoni"]<-NA
sortedSampleInfo$SpeciesCol[sortedSampleInfo$Species=="calliptera"]<-NA
sortedSampleInfo$SpeciesCol[sortedSampleInfo$Species=="alluaudi"]<-NA
sortedSampleInfo$SpeciesCol[sortedSampleInfo$Species=="gracilior"]<-"#000088"
sortedSampleInfo$SpeciesCol[sortedSampleInfo$Species=="pharyngalis"]<-"#1A1A88"
sortedSampleInfo$SpeciesCol[sortedSampleInfo$Species=="stappersi"]<-"#8B0000"
sortedSampleInfo$SpeciesCol[sortedSampleInfo$Species=="sp. Yaekama"]<-"#8B1A1A"

# Genus pattern
sortedSampleInfo$genusPattern<-as.factor(sortedSampleInfo$genusCol)
levels(sortedSampleInfo$genusPattern)<-rep(c(1:3),
                                           length.out=length(levels(sortedSampleInfo$genusPattern)))
sortedSampleInfo$genusPattern<-as.character(sortedSampleInfo$genusPattern)
sortedSampleInfo$genusPattern[is.na(sortedSampleInfo$genusPattern)]<-1


##### Prepare for branch colours #####

# Get lake colours for branches
clades<-list(which(sortedSampleInfo$lake_color=="darkblue"),
             which(sortedSampleInfo$lake_color=="#40A3DBB3"),
             which(sortedSampleInfo$lake_color=="blue"),
             which(sortedSampleInfo$lake_color=="deepskyblue"),
             which(sortedSampleInfo$lake_color=="red"),
             which(sortedSampleInfo$lake_color=="cornflowerblue"),
             which(sortedSampleInfo$lake_color=="#FF6304"),
             which(sortedSampleInfo$lake_color=="#FF9F00"),
             which(sortedSampleInfo$lake_color=="#F4C430"),
             which(sortedSampleInfo$lake_color=="black"&
                     sortedSampleInfo$Lake!="Malawi"&
                     sortedSampleInfo$Species!="calliptera"&
                     sortedSampleInfo$Species!="burtoni"),
             which(grepl("Nile",sortedSampleInfo$group)),
             which(sortedSampleInfo$Species=="stappersi"|
                     sortedSampleInfo$Species=="sp. Yaekama"))

treeInfo<-as_tibble(tree) %>% groupOTU(clades)

g1<-as(tree, 'phylo4')

d = data.frame("color"=treeInfo$group[nodeId(g1, "tip")])
rownames(d) = tree$tip.label

g2 <- phylo4d(g1, d)

rNodeData <- data.frame("color" = treeInfo$group[nodeId(g1, "internal")] ,
                        row.names = nodeId(g1, "internal"))



##### Genus colour palette ####
library(colorBlindness)
cols<-read.table("D:/Dropbox/R/colourPalette.txt",header=T)
pal <- rgb(cols$r, cols$g, cols$b)

pal<-pal[1:length(levels(factor(sortedSampleInfo$simplifiedEcology)))]

palGenera<-SteppedSequential5Steps[(1:length(SteppedSequential5Steps) %% 2) == 0]
palGenera<-c("black","darkgrey","grey","lightgrey","red","deeppink4","deeppink1","#B22C2C","#E57E7E","#99540FCC","#CC8E51","#FFD8B2",
             "#85B22C","#C3E57E","#0F6B99","#51A3CC","#B2E5FF","darkblue","#422CB2","#8F7EE5",
             "yellow","#FFFF0090","lightyellow")

# Remove genera outside of Victoria Radiation
sortedSampleInfo$genusCol[!(sortedSampleInfo$Lake=="Victoria"|sortedSampleInfo$Lake=="Kyoga")]<-NA



###### ggplot tree for Fig. 3 #######

p<-ggtree(g2,layout="circular",size=0.2,aes(col=factor(color))) %<+% sortedSampleInfo + 
  scale_colour_manual(name = "Lake",
                      values=c("white","darkblue","gold","blue","deepskyblue","red","cornflowerblue",
                               "#FF6304","#FF9F00","#F4C430","black","darkblue","darkred"),
                      labels=c("","Kivu","Sudan/Egypt/Turkana/Boukou","Edward",
                               "Saka/Mpanga River","Kagera","Albert","Kyoga","Victoria","Nabugabo/Ikimba/Burigi","other","Upper Nile","Congo"),
                      guide=guide_legend(ncol=1,keywidth = 0.5, keyheight = 0.5,order=1)) +
  new_scale_colour() +
  geom_tippoint(aes(color=SpeciesCol), size=0.4) +
  scale_colour_manual(
    values=levels(factor(sortedSampleInfo$SpeciesCol)),
    na.value = NA,
    labels = levels(factor(sortedSampleInfo$Species)),guide="none") +
  geom_nodepoint(size=ifelse(tree$node.label==100,0.5,ifelse(tree$node.label>80,0.25,0)), 
                 color="black", alpha=1) + 
  new_scale_fill() +
  new_scale_color() +
  geom_fruit(geom=geom_tile,width=0.005,
             mapping=aes(fill=genusCol),
             alpha=1,offset = 0.0,size = 0.07) +
  geom_tiplab2(aes(label=NA),col="lightgrey", align=T, 
               linetype = "12", linesize = 0.2, 
               size=1, offset=0.005, hjust=0)  +
  scale_fill_manual(
    values=alpha(palGenera,alpha = 0.4),
    labels=levels(factor(sortedSampleInfo$genusCol)),
    na.value = "white",
    guide=guide_legend(ncol=3,keywidth = 0.35, keyheight = 0.13,
                       order=2,title = "Genus")) +
  new_scale_fill() +
  new_scale_color() +
  geom_fruit(geom=geom_tile,width=0.002,
             mapping=aes(fill=lake_color),
             alpha=1,offset = 0.03,size = 0.04) +
  scale_fill_manual(values=levels(as.factor(sortedSampleInfo$lake_color)),
                    guide="none") +
  scale_color_manual(values=levels(as.factor(sortedSampleInfo$lake_color)),
                     guide="none") +
  theme(legend.justification = c(0.5,0.5), legend.position = c(0.52,0.55),
        plot.title = element_text(size = 6, face = "bold"),
        legend.title=element_text(size=4), 
        legend.text=element_text(size=3.2),
        legend.margin = margin(-0.3,0,0,0, unit="cm")) 

# Plot the tree rotated
rotate_tree(p,-148)



###### ggtree plot for Fig. S2 #######

p<-ggtree(g2,layout="circular",size=0.2,aes(col=factor(color))) %<+% sortedSampleInfo + 
  scale_colour_manual(name = "Lake",
                      values=c("white","darkblue","gold","blue","deepskyblue","red","cornflowerblue",
                               "#FF6304","#FF9F00","#F4C430","black","darkblue","darkred"),
                      labels=c("","Kivu","Sudan/Egypt/Turkana/Boukou","Edward",
                               "Saka/Mpanga River","Kagera","Albert","Kyoga","Victoria","Nabugabo/Ikimba/Burigi","other","Upper Nile","Congo"),
                      guide=guide_legend(ncol=1,keywidth = 0.5, keyheight = 0.5,order=1)) +
  new_scale_colour() +
  geom_tippoint(aes(color=SpeciesCol), size=0.4) +
  scale_colour_manual(
    values=levels(factor(sortedSampleInfo$SpeciesCol)),
    na.value = NA,
    labels = levels(factor(sortedSampleInfo$Species)),guide="none") +
  geom_nodepoint(size=ifelse(tree$node.label==100,0.5,ifelse(tree$node.label>80,0.25,0)), 
                 color="black", alpha=1) + 
  new_scale_fill() +
  geom_fruit(geom=geom_tile,width=0.005,
             mapping=aes(fill=genusCol),
             alpha=1,offset = 0.0,size = 0.07) +
  scale_fill_manual(
    values=colors(length(levels(factor(sortedSampleInfo$genusCol)))),
    labels=levels(factor(sortedSampleInfo$genusCol)),
    na.value = "white",
    guide=guide_legend(ncol=3,keywidth = 0.35, 
                       keyheight = 0.13,
                       order=2,title = "Genus")) +
  new_scale_fill() +
  geom_fruit(geom=geom_tile,width=0.005,
             mapping=aes(fill=simplifiedEcology),
             alpha=1,offset = 0.03,size = 0.07) +
  scale_fill_viridis_d(na.value="white",
                       guide=guide_legend(ncol=2,title="Ecomorph",
                                          keywidth = 0.35, keyheight = 0.13,order=3)) +
  new_scale_color() +
  new_scale_fill() +
  geom_fruit(geom=geom_tile,width=0.005,
             mapping=aes(fill=SpeciesCol),
             alpha=1,offset = 0.03,size = 0.07) +
  scale_fill_manual(values=levels(as.factor(sortedSampleInfo$SpeciesCol)),
                    guide="none") +
  scale_color_manual(values=levels(as.factor(sortedSampleInfo$SpeciesCol)),
                     guide="none") +
  new_scale_color() +
  geom_tiplab2(aes(label=NA,colour=SpeciesCol), align=T, 
               linetype = "dotted", linesize = 0.1, 
               size=1.3, offset=0.015, hjust=0)  +
  scale_color_manual(values=levels(as.factor(sortedSampleInfo$SpeciesCol)),
                     guide="none",na.value="lightgrey") +
  geom_tiplab2(aes(label=spLabel), align=T,color="black",
               size=1.7, offset=0.015, hjust=0,linetype=NA)  +
  theme(legend.justification = c(0.5,0.5), legend.position = c(0.52,0.55),
        plot.title = element_text(size = 7, face = "bold"),
        legend.title=element_text(size=6), 
        legend.text=element_text(size=5.2),
        legend.margin = margin(-0.3,0,0,0, unit="cm")) 

# Plot the tree rotated
rotate_tree(p,-148)




