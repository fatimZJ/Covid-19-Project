library("plotrix")

pdf("outputs/CCD_diagram.pdf", height = 7, width = 7)
plot(0, 0, pch = "", ylim = c(-2, 2), xlim = c(-2, 2), xlab = "", ylab = "", 
     axes = FALSE, asp = 1)
abline(v=0)
abline(h=0)

draw.circle(x = 0, y = 0, radius = 1, border = "blue")
verts <- 1/sqrt(2)
rect(-verts, -verts, verts, verts, border = "blue")

points(x = c(0, 0, 0, 1, -1, verts, -verts, verts, -verts), 
       y = c(0, 1, -1, 0 , 0, verts, verts, -verts, -verts),
       pch = "x", col = "red", cex = 2)
dev.off()

