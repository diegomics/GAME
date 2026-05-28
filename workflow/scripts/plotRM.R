#!/usr/bin/env Rscript
###############################################################################
# plotRM.R - Create repeat landscape plots from RepeatMasker data
#
# by Diego De Panis, 2025
# This script is part of the GAME pipeline
# note: AI tools may have been used to improve, clean and/or comment this version of the code
# 
# Usage: Rscript plotRM.R <configFile> <out.svg>
#
# Config file format: sample_id<TAB>path_to_tsv
###############################################################################

args <- commandArgs(trailingOnly=TRUE)

if (length(args) < 2) {
  stop("Error: 2 arguments needed!\nUsage: Rscript plotRM.R <configFile> <out.svg>", 
       call.=FALSE)
}

###############################################################################
# Repeat class definitions - ORDER IS IMPORTANT FOR STACKING
# This matches the order used by RepeatMasker's createRepeatLandscape.pl
###############################################################################

# Main landscape categories (bottom to top in stacked plot)
landscapeOrder <- c(
  'Unknown',
  'Other',
  'DNA.Academ',
  'DNA.CMC',
  'DNA.Crypton',
  'DNA.Ginger',
  'DNA.Harbinger',
  'DNA.hAT',
  'DNA.Kolobok',
  'DNA.Maverick',
  'DNA',
  'DNA.Merlin',
  'DNA.MULE',
  'DNA.P',
  'DNA.PiggyBac',
  'DNA.Sola',
  'DNA.TcMar',
  'DNA.Transib',
  'DNA.Zator',
  'DNA.Dada',
  'RC.Helitron',
  'LTR.DIRS',
  'LTR.Ngaro',
  'LTR.Pao',
  'LTR.Copia',
  'LTR.Gypsy',
  'LTR.ERVL',
  'LTR',
  'LTR.ERV1',
  'LTR.ERV',
  'LTR.ERVK',
  'LINE.L1',
  'LINE',
  'LINE.RTE',
  'LINE.CR1',
  'LINE.Rex.Babar',
  'LINE.L2',
  'LINE.Proto2',
  'LINE.LOA',
  'LINE.R1',
  'LINE.Jockey.I',
  'LINE.Dong.R4',
  'LINE.R2',
  'LINE.CRE',
  'PLE',
  'Retroposon.SVA',
  'SINE',
  'SINE.5S',
  'SINE.7SL',
  'SINE.Alu',
  'SINE.tRNA',
  'SINE.tRNA.Alu',
  'SINE.tRNA.RTE',
  'SINE.RTE',
  'SINE.Deu',
  'SINE.tRNA.V',
  'SINE.MIR',
  'SINE.U',
  'SINE.tRNA.7SL',
  'SINE.tRNA.CR1'
)

# Categories to EXCLUDE from landscape (used in pie chart only)
excludeCategories <- c(
  'Simple_repeat',
  'Satellite',
  'Structural_RNA',
  'Low_complexity'
)

# Color palette matching RepeatMasker canonical colors
customCols <- c(
  "Unknown"="#999999",
  "Other"="#4D4D4D",
  "DNA.Academ"="#FF0000",
  "DNA.CMC"="#FF200B",
  "DNA.Crypton"="#FF3115",
  "DNA.Ginger"="#FF3D1E",
  "DNA.Harbinger"="#FF4825",
  "DNA.hAT"="#FF512D",
  "DNA.Kolobok"="#FF5A34",
  "DNA.Maverick"="#FF623B",
  "DNA"="#FF6A42",
  "DNA.Merlin"="#FF7149",
  "DNA.MULE"="#FF7850",
  "DNA.P"="#FF7F57",
  "DNA.PiggyBac"="#FF865E",
  "DNA.Sola"="#FF8D65",
  "DNA.TcMar"="#FF936C",
  "DNA.Transib"="#FF9972",
  "DNA.Zator"="#FF9F79",
  "DNA.Dada"="#FFCFBC",
  "RC.Helitron"="#FF00FF",
  "LTR.DIRS"="#006400",
  "LTR.Ngaro"="#197214",
  "LTR.Pao"="#2A8024",
  "LTR.Copia"="#3A8F33",
  "LTR.Gypsy"="#489E42",
  "LTR.ERVL"="#57AE51",
  "LTR"="#65BD61",
  "LTR.ERV1"="#73CD70",
  "LTR.ERV"="#81DD80",
  "LTR.ERVK"="#90ED90",
  "LINE.L1"="#00008B",
  "LINE"="#251792",
  "LINE.RTE"="#38299A",
  "LINE.CR1"="#483AA2",
  "LINE.Rex.Babar"="#554BAA",
  "LINE.L2"="#625CB1",
  "LINE.Proto2"="#6E6DB9",
  "LINE.LOA"="#797EC0",
  "LINE.R1"="#848FC8",
  "LINE.Jockey.I"="#8FA1CF",
  "LINE.Dong.R4"="#99B3D7",
  "LINE.R2"="#A3C5DE",
  "LINE.CRE"="#C1D9FF",
  "PLE"="#ACD8E5",
  "Retroposon.SVA"="#FF4D4D",
  "SINE"="#9F1FF0",
  "SINE.5S"="#A637F1",
  "SINE.7SL"="#AD49F2",
  "SINE.Alu"="#B358F3",
  "SINE.tRNA"="#B966F4",
  "SINE.tRNA.Alu"="#BF74F4",
  "SINE.tRNA.RTE"="#C481F5",
  "SINE.RTE"="#C98EF6",
  "SINE.Deu"="#CE9BF7",
  "SINE.tRNA.V"="#D3A7F7",
  "SINE.MIR"="#D7B4F8",
  "SINE.U"="#DFCDF9",
  "SINE.tRNA.7SL"="#E2D9F9",
  "SINE.tRNA.CR1"="#E5E5F9"
)

###############################################################################
# Parse config file and create combined data table
###############################################################################
parseConfig <- function(configIN) {
  final_tab <- data.frame(
    Kimura=character(), 
    Per=numeric(), 
    Category=character(), 
    Sample=character(),
    stringsAsFactors=FALSE
  )
  
  # Loop over samples in config file
  for (i in seq(1, nrow(configIN), 1)) {
    sID <- configIN$V1[i]
    filePath <- configIN$V2[i]
    
    # Read TSV data
    curTab <- read.table(filePath, sep="\t", header=TRUE, check.names=FALSE)
    
    # Process each repeat class column (skip first column which is Divergence)
    for (j in seq(2, ncol(curTab), 1)) {
      category <- colnames(curTab)[j]
      
      # Convert slashes to dots to match color scheme
      category <- gsub("/", ".", category)
      
      df <- data.frame(
        Kimura = as.character(curTab[[1]]),  # Divergence column
        Per = as.numeric(curTab[[j]]),       # Percentage for this class
        Category = category,
        Sample = sID,
        stringsAsFactors = FALSE
      )
      
      final_tab <- rbind(final_tab, df)
    }
  }
  
  return(final_tab)
}

###############################################################################
# Main execution
###############################################################################

# Load required library
if (!require("ggplot2", quietly=TRUE)) {
  stop("Error: ggplot2 package is required but not installed", call.=FALSE)
}

# Read config file
confFile <- read.table(args[1], header=FALSE, sep="\t", stringsAsFactors=FALSE)

# Parse data
cat("Parsing repeat landscape data...\n")
repData <- parseConfig(confFile)

# Filter out non-landscape categories
repData <- repData[!repData$Category %in% excludeCategories, ]

# Set factor levels for proper stacking order (reverse for bottom-to-top)
# Only include categories that are actually present in the data
presentCategories <- intersect(rev(landscapeOrder), unique(repData$Category))
repData$Category <- factor(repData$Category, levels=presentCategories)

# Convert Kimura to numeric for proper x-axis ordering
repData$Kimura <- as.numeric(repData$Kimura)

# Remove any NA values
repData <- repData[complete.cases(repData), ]

# Report summary
cat(sprintf("Found %d categories across %d samples\n", 
            length(unique(repData$Category)), 
            length(unique(repData$Sample))))

# Create plot
cat("Creating SVG plot...\n")
svg(args[2], height=10, width=12)

p <- ggplot(repData, aes(x=Kimura, y=Per, fill=Category)) +
  geom_col(position="stack", width=1) + 
  scale_x_continuous(breaks=seq(0, 50, by=5), limits=c(0, 50)) +
  scale_y_continuous(expand=c(0, 0)) +
  scale_fill_manual(values=customCols, na.value="#CCCCCC", drop=TRUE) + 
  labs(
    title="Interspersed Repeat Landscape",
    x="Kimura substitution level (CpG adjusted)",
    y="Percent of genome"
  ) +
  facet_wrap(~Sample, ncol=3) + 
  theme_bw() +
  theme(
    legend.position="bottom",
    legend.key.size=unit(0.4, "cm"),
    legend.text=element_text(size=8),
    panel.grid.minor=element_blank(),
    plot.title=element_text(hjust=0.5, size=14, face="bold")
  ) + 
  guides(fill=guide_legend("", ncol=6, byrow=TRUE))

print(p)
dev.off()

cat("Done! Plot saved to:", args[2], "\n")
