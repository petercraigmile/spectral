
#tapers <- dpss.tapers(32, 2, NW=2)
#multitaper.effective.bw(tapers)


require(spectral)

#Fig 5 of Walden, McCoy and Percival (1995)

plot(1:3, sapply(1:3, function (k)
                 multitaper.effective.bw(dpss.tapers(32, k, NW=2))),
     xlab="", ylab="", ylim=c(0,0.25), xlim=c(1,7), type="l")
     
lines(1:5, sapply(1:5, function (k)
                 multitaper.effective.bw(dpss.tapers(32, k, NW=3))))

lines(1:7, sapply(1:7, function (k)
                 multitaper.effective.bw(dpss.tapers(32, k, NW=4))))
