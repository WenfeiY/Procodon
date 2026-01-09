#Creating syntenic plot
#@author: Sam_Y

library(circlize)

#  Load  data
genome_data <- read.csv("Ng.chr_df.txt", stringsAsFactors = FALSE)
genome_data$start <- as.numeric(genome_data$start)
genome_data$end <- as.numeric(genome_data$end)
links_data <- read.csv("Ng.sc_link.txt", stringsAsFactors = FALSE)
links_data$start1 <- as.numeric(links_data$start1)
links_data$end1 <- as.numeric(links_data$end1)
links_data$start2 <- as.numeric(links_data$start2)
links_data$end2 <- as.numeric(links_data$end2)

#  Draw figure
circos.clear()
circos.par(
  gap.degree = 5, 
  track.margin = c(0, 0),
  cell.padding = c(0, 0, 0, 0),
  start.degree = 90
)

circos.initialize(
  factors = genome_data$chr,
  xlim = cbind(genome_data$start, genome_data$end)
)

circos.track(
  ylim = c(0, 4),
  track.height = 0.1,
  bg.border = NA,

  panel.fun = function(x, y) {
    sector.index = get.cell.meta.data("sector.index")
    len = get.cell.meta.data("xlim")[2]
    
    species = ifelse(grepl("A", sector.index), "SpeciesA", "SpeciesB")
    col = ifelse(species == "SpeciesA", "#FF000080", "#0000FF80")
    
    circos.rect(0, 0, len, 1, col = col, border = col)
    
    label_pos = ifelse(species == "SpeciesA", 1.2, -0.2)
    circos.text(
      x = len/2,
      y = label_pos,
      labels = sector.index,
      facing = "downward",
      adj = c(0.5, 0.5),
      cex = 0.8
    )
  }
)

for(i in 1:nrow(links_data)) {
  circos.link(
    sector.index1 = links_data$chr1[i],
    point1 = c(links_data$start1[i], links_data$end1[i]),
    sector.index2 = links_data$chr2[i],
    point2 = c(links_data$start2[i], links_data$end2[i]),
    col = adjustcolor("purple", alpha.f = 0.6),
    border = NA,
    h.ratio = 0.5
  )
}

speciesA_chr = genome_data$chr[genome_data$species == "SpeciesA"]
speciesB_chr = genome_data$chr[genome_data$species == "SpeciesB"]

circos.text(
  x = mean(range(genome_data$end[1:3])),
  y = 0.5,
  labels = "Species A",
  sector.index = speciesA_chr[1],
  track.index = 1,
  cex = 1.2,
  col = "red",
  facing = "bending"
)

circos.text(
  x = mean(range(genome_data$end[4:5])),
  y = -0.5,
  labels = "Species B",
  sector.index = speciesB_chr[1],
  track.index = 1,
  cex = 1.2,
  col = "blue",
  facing = "bending"
)

legend("bottomright", 
       legend = c("Species A", "Species B", "Homologous Genes"),
       fill = c(adjustcolor("red", 0.5), adjustcolor("blue", 0.5), "purple"), 
       border = NA,
       cex = 0.9)

title("Genomic Synteny - Single Ring Visualization", cex.main = 1.2)
