####5###10########20########30########40########50########60########70########80
#
# This script creates figure 5, a panel of figs and fig wasps time series. To
# run this,
# 1. direct the file path from where you'll be working on line 16
# 2. Make sure the file path has two folders, one called "Data" (see line 17)
#    and the other called "Figures" (see line 18)
# 3. Make sure that "multifig_er.csv" is found in the Data folder
#
####5###10########20########30########40########50########60########70########80

# 0. Load libraries
	library(package = "viridisLite")

# 1. Set directories
	wd <- file.path("/Volumes/GoogleDrive/My Drive/Manuscripts/ObCoevo")
	dat_dir <- file.path(wd, "Data")
	fig_dir <- file.path(wd, "Figures")

# 2. Import data
	mp <- data.matrix(read.csv(file = file.path(dat_dir, "multifig_er.csv"), header = F))

# 3. Write functions

# 4. Prepare data and create varibles
	colnames(mp) <- c("step", "p1", "p2", "y", "z", "init", "eps")
	col_vec <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c")
	eps_vec <- unique(mp[,"eps"])
	eps_n <- length(eps_vec)
	init_vec <- unique(mp[,"init"])
	init_n <- length(init_vec)

# 5. Final figure 5
	png(filename = file.path(fig_dir, "Fig5_2pop_s6_r2_b02_multfigs2000.png"), res = 300, units = "in", width = 6, height = 7)
	par(mfrow = c(5, 3), mar = rep(x = 1/3, times = 4), oma = c(4, 4.5, 0.5, 0.5))
	for (i in 1:eps_n) {
			mp_i_ind <- which(mp[,"eps"] == eps_vec[i])
			mp_i <- mp[mp_i_ind,]
		for (j in 1:init_n) {
			plot(x = NA, type = "n", xlim = c(0, 2000), ylim = c(0, 1), ann = F, axes = F, frame = T)
			mp_j_ind <- which(mp_i[,"init"] == init_vec[j])
			mp_j <- mp_i[mp_j_ind,]
			lines(x = mp_j[,"step"], y = mp_j[,"p1"], col = col_vec[1], lwd = 1.5)
			lines(x = mp_j[,"step"], y = mp_j[,"p2"], col = col_vec[2], lwd = 1.5)
			lines(x = mp_j[,"step"], y = mp_j[,"y"], col = col_vec[3], lwd = 1.5)
			lines(x = mp_j[,"step"], y = mp_j[,"z"], col = col_vec[4], lwd = 1.5)
			if (j == 1) axis(side = 2, las = 1, at = seq(from = 0, to = 1, by = 0.25))
			if (i == eps_n) {
				axis(side = 1, at = pretty(seq(from = 0, to = 2000, by = 500)))
			}
			if (i == 1 & j == 1) {
				legend(x = "topright", col = col_vec, lty = 1, legend = c("A, wasp local adaptation locus", "B, wasp preference locus", "C, plant odd generations", "C, plant even generations"), bty = "n", cex = 0.75, lwd = 1.5)
			}
		}
	}
	mtext(side = 1, text = "Pollinator generations", outer = T, line = 2.5)
	mtext(side = 2, text = "Allele frequency in total populations", outer = T, line = 3)
	dev.off()