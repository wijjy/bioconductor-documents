## How to do animations in R the simple way

Using ordinary R functions and `imagemagick`.


```{r}
im <- function(i,prev) {
  jpeg(paste("a",i,".jpeg",sep=""))
  nv <- runif(2,-1,1)
  print(nv)
  here <- prev+nv
  print(here)
  ##bounce
  here[here<0] <- -here[here<0]
  here[here>20] <- 40-here[here>20]
  
  plot(here[1],here[2],xlim=c(0,20),ylim=c(0,20),xlab="",ylab="",axes=F)
  dev.off()
  here
}

here <- runif(2,0,10)
for (i in 1:100)  here <- im(i,here)

```{bash}
/usr/bin/convert  -page a1.jpeg -page a2.jpeg a3.jpeg a4.jpeg a5.jpeg a6.jpeg a7.jpeg a8.jpeg a9.jpeg a10.jpeg a.gif
```
