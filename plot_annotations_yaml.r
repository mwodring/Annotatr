 suppressPackageStartupMessages({
        library(GenomicRanges)
        library(Gviz)
        library(data.table)
        library(yaml)
    })

yaml_file <- "plot_anot.yaml"
yaml_data <- read_yaml(yaml_file)

plot_width <- yaml_data["width"][[1]]
plot_height <- yaml_data["height"][[1]]
contig_length <- yaml_data["contig_length"][[1]]
contig_name <- yaml_data$contig_name
plot_name <- yaml_data["plot_name"][[1]]
bg_file <-  yaml_data["bg_file"][[1]]

orf_csv <- yaml_data["orf_csv"]
orf_data <- read.csv(paste0(orf_csv))
anot_csv <- yaml_data["anot_csv"]
anot_data <- read.csv(paste0(anot_csv))
phob_csv <- yaml_data["phob_csv"]

if (identical(phob_csv, "")) {
    phob <- TRUE 
    phob_data <- read.csv(paste0(phob_csv))
    } else {
    phob <- FALSE 
    }

ORFs_starts <- orf_data$start
ORFs_ends <- orf_data$end
ORFs_widths <- ORFs_ends - ORFs_starts
ORFs_ids <- orf_data$name

aa_starts <- anot_data$Minimum
aa_ends <- anot_data$Maximum
aa_ids <- anot_data$Name

aa_ids <- replace(aa_ids, aa_ids=="", "N/A")

orf_shift_aa <- rep(ORFs_starts[1], times = length(aa_starts))
transform_aa <- rep(c(3), times = length(aa_starts))
aa_starts <- aa_starts * transform_aa
aa_ends <- aa_ends * transform_aa
aa_widths <- aa_ends - aa_starts

if (phob) {
phob_starts <- as.numeric(phob_data$start)
phob_ends <- as.numeric(phob_data$end)
phob_ids <- character()

for (phob_type in phob_data$type) {
    if (identical(phob_type, "NC")) { 
    phob_ids <- append(phob_ids, "lawngreen", after = length(phob_ids))
    } else if (identical(phob_type, "TMR")) {
        "yes"
    phob_ids <- append(phob_ids, "magenta", after = length(phob_ids))
    } else if (identical(phob_type, "C")) {
    phob_ids <- append(phob_ids, "gold1", after = length(phob_ids))
    } else {
    print("Please use only NC, TMR and C in your .csv file.")
}
}

transform_phob <- rep(c(3), times = length(phob_starts))
phob_starts <- phob_starts * transform_phob
phob_ends <- phob_ends * transform_phob
phob_widths <- phob_ends - phob_starts
phobs <- GRanges(seqnames = Rle(contig_name, c(length(phob_ids))), 
                ranges = IRanges(phob_starts, width = phob_widths, names = phob_ids))
    }

whole_contig <- GRanges(seqnames = Rle(c(contig_name), c(1)), 
                ranges = IRanges(0, width = contig_length))
whole_aa <- GRanges(seqnames = Rle(c(contig_name), c(1)), 
                ranges = IRanges(ORFs_starts, width = ORFs_widths))

ORFs <- GRanges(seqnames = Rle(c(contig_name), c(length(ORFs_ids))), 
                ranges = IRanges(ORFs_starts, width = ORFs_widths, names = ORFs_ids))
prots <- GRanges(seqnames = Rle(c(contig_name), c(length(aa_ids))), 
                ranges = IRanges(aa_starts, width = aa_widths, names = aa_ids))

bedgraph_dt <- fread(bg_file, col.names = c('chromosome', 'start', 'end', 'value'))

options(ucscChromosomeNames=FALSE)

atrack <- AnnotationTrack(range = ORFs, 
                          name = plot_name, 
                          id = ORFs_ids)
atrack3 <- AnnotationTrack(range = prots,
                           name = "Domains",
                           id = aa_ids)

gtrack2 <- GenomeAxisTrack(whole_aa, littleTicks = TRUE, cex = 1, name = "AAs")

datrack <- DataTrack(range = bedgraph_dt, genome = whole_contig,
                     chromosome = contig_name,
                     name = "Coverage")
datrack2 <- DataTrack(range = bedgraph_dt, genome = whole_contig,
                      chromosome = contig_name,
                      name = "Line") 

displayPars(atrack) <- list(fontcolor.item = "black", 
                            col = "darkblue", fill = c("lightblue", "lightblue", "lightblue", "coral", "lightblue"),
                            showFeatureId = TRUE, background.title = "darkgray",
                            stacking = "squish")
displayPars(atrack3) <- list(fontcolor.item = "black", 
                            col = "midnightblue", fill = c("seagreen1", "turquoise3", "springgreen3"), detailsBorder.col = "blue",
                            showFeatureId = TRUE, background.title = "darkgray",
                            stacking = "dense", rotation.item = 90,
                            alpha = 0.75, alpha.title = 1)

displayPars(datrack) <- list(type = "gradient", 
                             gradient = c("mintcream", "lightskyblue1", "paleturquoise3", "lightsalmon", "orange", "orangered1"),
                             background.title = "darkgray", cex.axis = 1)
displayPars(datrack2) <- list(type = "a", alpha.title = 0, col= "black")

otrack <- OverlayTrack(trackList=list(datrack, datrack2), 
                                       name="Coverage", background.title = "darkgray")

if (phob) {
    atrack2 <- AnnotationTrack(range = phobs, name = "Membrane", id = phob_ids)
    displayPars(atrack2) <- list(fontcolor.item = "black", 
                            col = "darkblue", fill = phob_ids,
                            background.title = "darkgray",
                            stacking = "dense")
    }

if (phob) {
    toPlot <- list(atrack, otrack, gtrack2, atrack3, atrack2) 
    plot_sizes <- c(2,5,1,2,2)
    } else {
    toPlot <- list(atrack, otrack, gtrack2, atrack3)
    plot_sizes <- c(2,5,1,2)
    }

jpeg(paste0(plot_name, ".jpeg"), width = plot_width, height = plot_height)
                plotTracks(toPlot, add53= TRUE, stackHeight = 0.9, add = TRUE, sizes = plot_sizes)
                dev.off()
