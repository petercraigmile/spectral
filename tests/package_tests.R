

library(spectral)

y <- arima.sim(n=1024, model=list(ar=ar2.example), sd=0.2)

sp <- pgram(y)
plot(sp)


sp <- direct.spectral(y, "cosine", p=0.2)
plot(sp)




sp <- direct.spectral(y, "dpss", n.tapers=7)

plot(mean(sp))

plot(sp, 1)
