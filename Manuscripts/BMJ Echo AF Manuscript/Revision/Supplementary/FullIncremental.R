rm(list=ls())
setwd("X:/BMJ Echo AF Manuscript/Revision/Supplementary")

Data=data.frame(
cost=c(
2459,
4712,
2815,
5405,
1527,
2467,
1974,
3106,
2449,
4614,
2779,
5315,
1510,
2393,
1955,
3039,
1487,
2321,
1942,
2946
),

qaly=c(
13.60,
13.51,
14.27,
14.19,
9.12,
9.13,
9.94,
9.97,
13.61,
13.54,
14.27,
14.22,
9.12,
9.15,
9.95,
9.99,
9.13,
9.18,
9.95,
10.01
),
opt=c(
"wm50n",
"wm50t",
"wf50n",
"wf50t",
"wm65n",
"wm65t",
"wf65n",
"wf65t",

"rm50n",
"rm50t",
"rf50n",
"rf50t",
"rm65n",
"rm65t",
"rf65n",
"rf65t",

"dm65n",
"dm65t",
"df65n",
"df65t"
)
)

Data.m50 <- Data[c(1,2,9,10),]
Data.f50 <- Data[c(3,4,11,12),]
Data.m65 <- Data[c(5,6, 13,14, 17, 18),]
Data.f65 <- Data[c(7,8, 15,16, 19,20),]


jpeg("Fullincm50.jpg", height=500, width=500)

# males, age 50 
plot(
  Data.m50$qaly, Data.m50$cost,
  type="p", 
  pch=c(0, 15, 1, 16),
  xlab="Mean Health Benefit (QALYs)",
  ylab="Mean Cost (£)"
)

legend(
  "topright", 
  pch=c(0,15, 1,16),  
  legend=c("(X) Warf. No TTE",
           "(X) Warf. TTE",
           "(1) Riv. No TTE", 
           "(X) Riv. TTE" 
       )
)

dev.off()

jpeg("Fullincf50.jpg", height=500, width=500)


# females, age 50 
plot(
  Data.f50$qaly, Data.f50$cost,
  type="p", 
  pch=c(0, 15, 1, 16),
  xlab="Mean Health Benefit (QALYs)",
  ylab="Mean Cost (£)"
)

legend(
  "topright", 
  pch=c(0,15, 1,16),  
  legend=c("(X) Warf. No TTE",
           "(X) Warf. TTE",
           "(1) Riv. No TTE", 
           "(X) Riv. TTE" 
  )
)

dev.off()

jpeg("Fullincm65.jpg", height=500, width=500)


# males, age 65 
plot(
  Data.m65$qaly, Data.m65$cost,
  type="p", 
  pch=c(0, 15, 1, 16, 2, 17),
  xlab="Mean Health Benefit (QALYs)",
  ylab="Mean Cost (£)"
)

lines(Data.m65$qaly[c(5,6)], Data.m65$cost[c(5,6)], lwd=2, lty="dashed")

legend(
  "bottomright", 
  pch=c(0,15, 1,16, 2, 17, NA),
  lwd=c(rep(NA, 6), 2),
  lty=c(rep(NA, 6), "dashed"),
  legend=c("(X) Warf. No TTE",
           "(X) Warf. TTE",
           "(X) Riv. No TTE", 
           "(X) Riv. TTE",
           "(1) Dab. No TTE",
           "(2) Dab. TTE",
           "Frontier"
  )
)
dev.off()


# Females, age 65
jpeg("Fullincf65.jpg", height=500, width=500)

plot(
  Data.f65$qaly, Data.f65$cost,
  type="p", 
  pch=c(0, 15, 1, 16, 2, 17),
  xlab="Mean Health Benefit (QALYs)",
  ylab="Mean Cost (£)"
)

lines(Data.f65$qaly[c(5,6)], Data.f65$cost[c(5,6)], lwd=2, lty="dashed")

legend(
  "bottomright", 
  pch=c(0,15, 1,16, 2, 17, NA),
  lwd=c(rep(NA, 6), 2),
  lty=c(rep(NA, 6), "dashed"),
  legend=c("(X) Warf. No TTE",
           "(X) Warf. TTE",
           "(X) Riv. No TTE", 
           "(X) Riv. TTE",
           "(1) Dab. No TTE",
           "(2) Dab. TTE",
           "Frontier"
  )
)

dev.off()


