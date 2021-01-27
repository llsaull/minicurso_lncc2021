# Faz simulações do modelo usado em Zhang 2001
#  usando o Euler-Maruyama.

source("EulerMaruyama.R")

# Dados do artigo Zhang
lambda <- c(6.04, 8.90)
sigma  <- c(0.44, 0.63) 
r      <- c(1.5, -1.61)
mu     <- r + sigma^2/2
X0 = 70.5
p0 = c(0.5,0.5)

# Fazendo os gráficos dos caminhos
L  <- 2**12
T  <- 1
Dt <- T/L

# Fazendo a semente ser 15 para dar sempre o mesmo caminho
set.seed(15)

# Simulação
time <- seq(0,T,Dt)
res  <- runEulerMaruyama(X0,p0,mu,sigma,lambda,Dt,L)
X    <- res$X
A    <- res$A

par(mfrow=c(1,1))
plot(time,X,main="",
     ylab="",xlab="Time",type="l",col="red",lwd=2)
par(new=TRUE)
plot(time,A,ylab="",xlab="Time",type="s",lty=2,col="blue",lwd=2,axes=FALSE)
axis(4, at=c(1,2),las=1)

# Várias simulações:
S <- 30 
plot(time,X,main="",
     ylab="",xlab="Time",type="l",col="red",lwd=2)

# Contabilizando o número de caminhos que passam de 100 antes de cruzar 50.
cGoodPath <- 0
MIN       <- 50
MAX       <- 100
goodPath  <- function(X){
    if(min(X) > MIN && max(X) > MAX) return(TRUE)
    else if(max(X) > MAX && which.min(X) > which.max(X)) return(TRUE)
    else return(FALSE)
}
if(goodPath(X)){cGoodPath <- cGoodPath + 1 }

for( i in 1:S )
{
    X <- runEulerMaruyama(X0,p0,mu,sigma,lambda,Dt,L)$X
    lines(time,X,main="",
         ylab="",xlab="Time",type="l",col=rgb(0.5,0.1,0.5,0.2),lwd=2)
    
    if(goodPath(X)){cGoodPath <- cGoodPath + 1 }
} 
cat("prob:",cGoodPath/S)
