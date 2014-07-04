X1 <- rnorm(1000, 2, 2)
X2 <- rnorm(1000, 4, 1)
X3 <- rnorm(1000, 6)

Data.unsorted <- data.frame(x1=X1, x2=X2, x3=X3)

Data.sorted <- data.frame(x1=sort(X1), x2=sort(X2), x3=sort(X3))

plot(x2 - x1 ~ 1, data=Data.sorted, main ="Differences between x1 and x2", ylim=c(min(Data.sorted), max(Data.sorted)))
points(Data.unsorted[,2] - Data.unsorted[,1], col="grey")
abline(h=0, lwd=2, lty="dashed")
legend("topright", pch=c(1,1), col=c("black", "grey"), legend=c("Sorted", "Unsorted"))

# systematise this:

sd.ratios <- 10^(seq(-2, 2, by=0.01))
N.sd.ratios <- length(sd.ratios)
proportion.below.unsorted <- vector("numeric", N.sd.ratios)
proportion.below.sorted <- vector("numeric", N.sd.ratios)
var.diff.unsorted <- vector("numeric", N.sd.ratios)
var.diff.sorted <- vector("numeric", N.sd.ratios)

for (i in 1:N.sd.ratios){
    X1 <- rnorm(100000, 1, 1)
    X2 <- rnorm(100000, 2, 1 * sd.ratios[i])
    
    proportion.below.unsorted[i] <- length(which(X1 > X2))  / 100000
    proportion.below.sorted[i] <- length(which(sort(X1) > sort(X2))) / 100000
    var.diff.unsorted[i] <- var(X2-X1)
    var.diff.sorted[i] <- var(sort(X2) - sort(X1))
}

Data <- data.frame(unsorted= proportion.below.unsorted, sorted=proportion.below.sorted,  var.diff.srt= var.diff.sorted, var.diff.unsrt=var.diff.unsorted, ratio=sd.ratios)
plot(unsorted ~ ratio, data=Data, ylim=c(0,0.5), col="grey", type="l", log="x")
lines(sorted ~ ratio, data=Data, col="black")

plot(var.diff.srt ~ ratio, data=Data, type="l", col="black", lwd=2, log="xy", ylab="variance of differences", xlab="ratio of SDs")

lines(var.diff.unsrt ~ ratio, data=Data, lwd=2, col="black", lty="dashed")
legend("topleft", legend=c("sorted", "unsorted"), lwd=c(2,2), col=c("black", "black"), lty=c("solid", "dashed"))

save(Data, file="D:/Work/EchoAF/SimData.rData")
