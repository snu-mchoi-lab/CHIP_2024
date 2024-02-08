
power_n<-function(size){
    p0 = 1/40 #for calculate the critical region
    c<-qbinom(p0,size,0.025)
    p<-seq(0.01,0.4,0.01)
    pow<-1-dbinom(c,size,p)+dbinom(size-c,size,p)
    return(pow)
}

p<-seq(0.01,0.4,0.01)
print(p)
plot(p,power_n(235),type="l",col="red")
lines(p,power_n(300),lty=2,col="blue")
lines(p,power_n(400),lty=3, col="black")

text(0.23,0.93,"Red, sample size =235")
text(0.23,0.94,"Blue, sample size =300")
text(0.24,0.95,"Black, sample size =400")