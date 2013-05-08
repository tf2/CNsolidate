`vcbar` <- function(zlim, colV) 
 {
  n <- length(colV)
  plot(NA, xlim = c(0,1), ylim = c(1,n), axes = FALSE, xlab = "", ylab = "",
       main = "", xaxs = "i", yaxs = "i")
  abline(h = 1:n, col = colV, lwd = 1)
  axis(4, at = round(seq(1,n,length.out=5)), las = 2, cex.axis = 0.6,
       label = formatC(seq(zlim[1], zlim[2], length.out = 5),
         format = "f", digits = 2))
  box()
  mtext("Log-Ratio", side = 3, line = 0.5, cex = 0.7, font = 2,
        at = 0.9)
}

gllegend <-
function(Class.names)  {
if (length(Class.names) > 2)
{
  plot(NA, xlab = "", ylab = "", main = "", axes = FALSE,
       xlim = c(0,1), ylim = c(0,1), xaxs = "i", yaxs = "i")
  legend(x = "center", col = c("red","orange", "green","green4"), pch = 15, cex = 0.9,
         pt.cex = 1.2, bg = "white", legend = Class.names, bty = "n")
} else
{
 plot(NA, xlab = "", ylab = "", main = "", axes = FALSE,
       xlim = c(0,1), ylim = c(0,1), xaxs = "i", yaxs = "i")
  legend(x = "center", col = c("orange", "green"), pch = 15, cex = 0.8,
         pt.cex = 1.2, bg = "white", legend = Class.names, bty = "n")

}
}

