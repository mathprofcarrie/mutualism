####5###10########20########30########40########50########60########70########80
#
# This script creates figure 2, a time secies with host switching. To run this,
# 1. direct the file path from where you'll be working on line 14
# 2. Make sure the file path has two folders, one called "Data" (see line 15)
#    and the other called "Figures" (see line 16)
# 3. Make sure that "greya_switch.csv" is found in the Data folder
#
####5###10########20########30########40########50########60########70########80

# 0. Load libraries

# 1. Set directories
	wd <- file.path("/Volumes/GoogleDrive/My Drive/Manuscripts/ObCoevo")
	dat_dir <- file.path(wd, "Data")
	fig_dir <- file.path(wd, "Figures")

# 2. Import data
	gs <- data.matrix(read.csv(file = file.path(dat_dir, "greya_switch.csv"), header = F))

# 3. Write functions

# 4. Prepare data and create varibles
	colnames(gs) <- c("gen", "A", "B", "C")
	col_vec <- c("#a6cee3", "#1f78b4", "#b2df8a")

# 5. Plot figure
	png(file = file.path(fig_dir, "Fig2.png"), height = 5, width = 5, units = "in", res = 300)
		par(mar = c(4.5, 4.5, 0.1, 0.1))
		plot(x = gs[,"gen"], y = gs[,"A"], type = "l", lwd = 2, col = col_vec[1], xlab = "Generation", ylab = "Allele frequency in total population", las = 1, ylim = c(0, 1))
			lines(x = gs[,"gen"], y = gs[,"B"], lwd = 2, col = col_vec[2])
			lines(x = gs[,"gen"], y = gs[,"C"], lwd = "2", col = col_vec[3])
			legend(x = "topright", legend = c("A, moth local adaptation locus", "B, moth preference locus", "C, plant"), lwd = 2, col = col_vec, bty = "n")
	dev.off()