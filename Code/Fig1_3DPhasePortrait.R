####5###10########20########30########40########50########60########70########80
#
# This script creates figure 1, the 3-dimensional phase portrait. To run this,
# 1. direct the file path from where you'll be working on line 16
# 2. Make sure the file path has two folders, one called "Data" (see line 17)
#    and the other called "Figures" (see line 18)
# 3. Make sure that "ngreya.csv" is found in the Data folder
#
####5###10########20########30########40########50########60########70########80

# 0. Load libraries
	library(package = "plot3D")
	library(package = "viridisLite")

# 1. Set directories
	wd <- file.path("/Volumes/GoogleDrive/My Drive/Manuscripts/ObCoevo")
	dat_dir <- file.path(wd, "Data")
	fig_dir <- file.path(wd, "Figures")

# 2. Import data
	ng <- data.matrix(read.csv(file = file.path(dat_dir, "ngreya.csv"), header = F))

# 3. Write functions
	panelfirst <- function(pmat) {
		zmin <- 0
		XY <- trans3D(mat[rep_i, 3], mat[rep_i, 4], z = zmin, pmat = pmat)
		lines2D(XY$x, XY$y, pch = 16, col = "#00000022", cex = 0.5, add = TRUE, colkey = FALSE)
		xmin <- 0
		XY <- trans3D(x = xmin, y = mat[rep_i, 4], z = mat[rep_i, 5], pmat = pmat)
		lines2D(XY$x, XY$y, pch = 16, col = "#00000033", cex = 0.5, add = TRUE, colkey = FALSE)
		ymax <- 1
		XY <- trans3D(x = mat[rep_i, 3], y = ymax, z = mat[rep_i, 5], pmat = pmat)
		lines2D(XY$x, XY$y, pch = 16, col = "#00000033", cex = 0.5, add = TRUE, colkey = FALSE)
	}

# 4. Reformat data to long form
	replicates <- 100
	n_steps <- nrow(ng)
	step <- rep(x = ng[,1], times = replicates)
	thirds <- seq(from = 2, to = ncol(ng), by = 3)
	p1 <- as.vector(ng[, thirds])
	p2 <- as.vector(ng[, thirds + 1])
	q <- as.vector(ng[, thirds + 2])
	replicate <- rep(x = 1:replicates, each = n_steps)
	mat <- matrix(data = c(replicate, step, p1, p2, q), ncol = 5)
		colnames(mat) <- c("rep", "step", "p1", "p2", "q")

# 5. Final figure 1
	# Treajectory colored by time
	png(filename = file.path(fig_dir, "Fig1_col.png"), res = 300, units = "in", width = 7, height = 7)
		par(mfrow = c(1, 1), mar = c(6, 4, 3, 5))
		scatter3D(x = 0.5, y = 0.5, z = 0.5, xlab = "\n\nFrequency of allele A\nin insect population", ylab = "\n\nFrequency of allele B\n in insect population", zlab = "", xlim = c(0, 1), ylim = c(0, 1), zlim = c(0, 1), theta = 330, phi = 15, ticktype = "detailed", r = 50, border = NA, col = "white", bty = "f")
		mtext(side = 2, text = "\n\nFrequency of allele C\n in plant population", line = 1.75)
		for (i in 1:replicates) {
		rep_i <- which(mat[,1] == i)
		lines3D(x = mat[rep_i, 3], y = mat[rep_i, 4], z = mat[rep_i, 5], col = viridis(n = n_steps, alpha = 0.5), add = T, pch = 19, lwd = 1.5, colvar = log(c(1:n_steps)), colkey = list(plot = FALSE))
		}
		scatter3D(x = c(0, 1), y = c(0, 1), z = c(0, 1), add = T, pch = 19, col = "red")
		colkey(col = rev(viridis(n = n_steps)), clim = rev(c(1, 500)), add = T, clab = "Generation", length = 0.75, clog = TRUE)
	dev.off()