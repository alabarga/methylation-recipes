
library("devtools")

source_url("http://raw.github.com/nielshanson/mp_tutorial/master/downstream_analysis_r/code/venn_diagram2.r")
source_url("http://raw.github.com/nielshanson/mp_tutorial/master/downstream_analysis_r/code/venn_diagram3.r")
source_url("http://raw.github.com/nielshanson/mp_tutorial/master/downstream_analysis_r/code/venn_diagram4.r")


lista <- sample(rownames(EPIC.manifest.hg38), 3000)
list1 <- sample(lista, 1000)
list2 <- sample(lista, 1000)
list3 <- sample(lista, 1000)

plot.new()
v2 <- venn_diagram2(list1, list2, "lista I", "lista II")
plot.new()
v3<- venn_diagram3(list1, list2, list3, "lista I", "lista II", "lista III")