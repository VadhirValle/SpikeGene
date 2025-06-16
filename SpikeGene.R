library(seqinr)
library(dplyr)
library(ggplot2)

trad = c(
  UUU="F", UUC="F", UUA="L", UUG="L",
  UCU="S", UCC="S", UCA="S", UCG="S",
  UAU="Y", UAC="Y", UAA="STOP", UAG="STOP",
  UGU="C", UGC="C", UGA="STOP", UGG="W",
  CUU="L", CUC="L", CUA="L", CUG="L",
  CCU="P", CCC="P", CCA="P", CCG="P",
  CAU="H", CAC="H", CAA="Q", CAG="Q",
  CGU="R", CGC="R", CGA="R", CGG="R",
  AUU="I", AUC="I", AUA="I", AUG="M",
  ACU="T", ACC="T", ACA="T", ACG="T",
  AAU="N", AAC="N", AAA="K", AAG="K",
  AGU="S", AGC="S", AGA="R", AGG="R",
  GUU="V", GUC="V", GUA="V", GUG="V",
  GCU="A", GCC="A", GCA="A", GCG="A",
  GAU="D", GAC="D", GAA="E", GAG="E",
  GGU="G", GGC="G", GGA="G", GGG="G"
)

datos = data.frame(
  mutacion = character(),
  cambioCodon = character(),
  cambioAmino = character(),
  pos = integer(),
  gen = character()
)

file = read.fasta("sequence.txt", forceDNAtolower = FALSE)
file2 = read.fasta("582MEXB11519.fasta", forceDNAtolower = FALSE)
cat (length(file)%/%12, "vs", length(file2)/12, "secuencias \n")

nMut=1
for (i in seq_along(file)){
  if (i==2) next
  gen = file[[i]]
  info = attr(gen,"Annot")
  info = unlist(strsplit(info,"\\[|\\]|:|=|\\."));
  gene = info[which(info=="gene")+1]
  
  if (gene != "S") next
  
  cat ("Gen",i,gene,"\n")
  gen[gen=="T"] = "U"
  cat("Total de nucleótidos (Wuhan):", length(gen), "\n")
  
  for (j in seq(i, length(file2), 12)){
    gen2 = file2[[j]]
    gen2[gen2=="T"] = "U"
    
    if (length(gen) == length(gen2)){
      diff = which(gen != gen2)
      if (length(diff) > 0){ 
        cat("Se encontraron mutaciones en las posiciones:", diff, "\n")
        prevMut = ""
        for (pos in diff){
          ini = pos - (pos-1)%%3
          mutacion = paste(gen[pos], "to", gen2[pos], sep="")
          codOri = paste(gen[ini],gen[ini+1],gen[ini+2],sep="")
          codMut = paste(gen2[ini],gen2[ini+1],gen2[ini+2],sep="")
          codonChange = paste(codOri,"to",codMut,sep="")
          nCod = ((pos-1)%/%3) + 1
          aminoChange = paste(trad[codOri],nCod,trad[codMut],sep="")
          
          if (!is.na(trad[codMut]) && trad[codOri]!=trad[codMut] && prevMut != aminoChange){
            cat(mutacion, codonChange, aminoChange, nCod, gene, "\n")
            datos[nMut, ] = list(mutacion, codonChange, aminoChange, nCod, gene)
            nMut = nMut + 1
          }
          prevMut = aminoChange
        }
      }
    }
  }
}
n_sec <- length(file2) / 12

# GRAFICA A: MUTACIONES ≥50%

frecuencia_mut <- datos %>%
  count(mutacion, name = "freq") %>%
  filter(freq / n_sec >= 0.50)

frecuencia_amino <- datos %>%
  count(cambioAmino, name = "freq") %>%
  filter(freq / n_sec >= 0.50)

datos_filtrados <- datos %>%
  filter(mutacion %in% frecuencia_mut$mutacion,
         cambioAmino %in% frecuencia_amino$cambioAmino)

ggplot(datos_filtrados, aes(x = mutacion)) +
  geom_bar(fill = "indianred") +
  labs(title = "A: Mutaciones de nucleótido (≥50% de las secuencias)",
       x = "Mutación",
       y = "Frecuencia") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

# GRAFICA B: CAMBIO AMINO ≥50

dfgraph <- datos %>%
  group_by(cambioAmino) %>%
  summarise(
    mutacion = first(mutacion),
    cambioCodon = first(cambioCodon),
    pos = first(pos),
    gen = first(gen),
    cuenta = n()
  ) %>%
  filter(cuenta >= 50)

ggplot(dfgraph) +
  aes(x = cambioAmino, y = cuenta, fill = cambioAmino, label = cuenta) +
  ggtitle("B: Cambios de aminoácidos mejorada (≥50 ocurrencias)") +
  labs(x = "Cambio", y = "Frecuencia", fill = "Frecuencia") +
  geom_bar(stat = "identity") +
  geom_text(stat = "identity", vjust = 1.5) +
  facet_grid(~gen, scales = "free", space = "free_x") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# ANALISIS: Mutaciones de Interés (VoC / VoI)

mutaciones_voc <- data.frame(
  cambioAmino = c("D614G", "N501Y", "E484K", "T478K", "P681R", "L452R", "K417N", "K417T"),
  variante_asociada = c(
    "Alpha/Delta/Ómicron", "Alpha/Beta/Gamma/Ómicron", "Beta/Gamma", "Delta/Ómicron", 
    "Delta", "Delta", "Beta/Ómicron", "Gamma"),
  stringsAsFactors = FALSE
)

datos_convoc <- datos %>%
  inner_join(mutaciones_voc, by = "cambioAmino")

tabla_resumen_voc <- datos_convoc %>%
  count(cambioAmino, variante_asociada, name = "ocurrencias") %>%
  arrange(desc(ocurrencias))

cat("Mutaciones de interés encontradas:\n")
print(tabla_resumen_voc)

# Paso 4: Gráfica de frecuencia de VoC (variants of concern)
ggplot(datos_convoc, aes(x = cambioAmino, fill = variante_asociada)) +
  geom_bar() +
  labs(title = "C: Mutaciones VoC en variantes mexicanas",
       x = "Cambio de aminoácido",
       y = "Frecuencia",
       fill = "Variante asociada") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))