####5###10########20########30########40########50########60########70########80
#
# This script creates figure 3, a phase portrait of p_2 and y_1. To run this,
# 1. direct the file path from where you'll be working on line 16
# 2. Make sure the file path has two folders, one called "Data" (see line 17)
#    and the other called "Figures" (see line 18)
# 3. Make sure that "multifig_er.csv" is found in the Data folder
#
####5###10########20########30########40########50########60########70########80

# 0. Load libraries
	library(package = "deSolve")
	library(package = "phaseR")

# 1. Set directories
	wd <- file.path("/Volumes/GoogleDrive/My Drive/Manuscripts/ObCoevo")
	fig_dir <- file.path(wd, "Figures")

# 2. Write diff. eq. function compatable with deSolve and phaseR
	eqn <- function(t, y, parms) {
		with(as.list(c(y, parms)), {
			q2 <- 1 - p2
			y2 <- 1 - y1
			dp2 <- 3*epsillon*p2*q2*(y1 - y2)
			dy1 <- 2*beta*epsillon*y1*y2*(p2 - q2)
			return(list(c(dp2, dy1)))
		})
	}

# 3. Set up parameter and evaluation conditions
	parms <- c(epsillon = 0.1, beta = 0.2)
	lims <- c(0, 1)

# 4. Find equilibria points for the phase plane
	y_star_interior <- findEquilibrium(deriv = eqn, y0 = c(0.5, 0.5), parameters = parms, state.names = c("p2", "y1"))$ystar
	y_star_trivial <- findEquilibrium(deriv = eqn, y0 = c(0.1, 0.1), parameters = parms, state.names = c("p2", "y1"))$ystar
	y_star_fixed <- findEquilibrium(deriv = eqn, y0 = c(0.9, 0.9), parameters = parms, state.names = c("p2", "y1"))$ystar

# 5. Create figure
	pdf(file = file.path(fig_dir, "PhasePlane.pdf"), height = 5, width = 5, units = "in")

	par(xpd = T, mar = c(3.5, 3.5, 1, 1)) # xpd = T is b/c the eqilibirum points are at the flush asix limits
	plot(x = NA, xlim = lims, ylim = lims, xaxs = "i", yaxs = "i", las = 1, ann = F)
		mtext(side = 1, text = bquote(p[2]), line = 2.25)
		mtext(side = 2, text = bquote(y[1]), line = 2.5, las = 2)
		grid_intervals <- (1:9)/10
		segments(x0 = grid_intervals, x1 = grid_intervals, y0 = rep(0, times = 9), y1 = rep(1, times = 9), col = "grey75")
		segments(x0 = rep(0, times = 9), x1 = rep(1, times = 9), y0 = grid_intervals, y1 = grid_intervals , col = "grey75")
	
	# Dumping the output of phaseR's functions 
		trash_out <- flowField(deriv = eqn, xlim = lims, ylim = lims, parameters = parms, state.names = c("p2", "y1"), add = T, arrow.type = "proportional", arrow.head = 0.025, col = "darkgreen", points = 21, frac = 4/4)
		trash_out <- drawManifolds(deriv = eqn, y0 = c(0.5, 0.1), xlim = lims, ylim = lims, parameters = parms, state.names = c("p2", "y1"), add.legend = FALSE, col = c("#845DC2", "#5594B9"), lwd = 1.5)
	rm(trash_out)
	
		points(x = y_star_interior[1,], y = y_star_interior[2,], col = "red", bg = "white", lwd = 1.5, cex = 1.25, pch = 21)
		points(x = y_star_trivial[1,], y = y_star_trivial[2,], col = "red", bg = "red", lwd = 1.5, cex = 1.25, pch = 21)
		points(x = y_star_fixed[1,], y = y_star_fixed[2,], col = "red", bg = "red", lwd = 1.5, cex = 1.25, pch = 21)
	dev.off()