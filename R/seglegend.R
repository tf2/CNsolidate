seglegend <-
function(Method)  {
  plot(NA, xlab = "", ylab = "", main = "", axes = FALSE,
       xlim = c(0,1), ylim = c(0,1), xaxs = "i", yaxs = "i")
  legend(x = "center", col = c("grey","black"), pch=c(19,1), lty=1,lwd=c(0,1.5), cex = 0.9,
         pt.lwd = c(0.9,0), bg = "white", legend = c("data",Method), bty = "n")
}
