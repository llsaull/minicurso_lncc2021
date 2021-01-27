# Este script tenta resolver o problema 
#  apresentado por Zhang 2001 usando o 
#  método da cadeia de Markov aproximada.

# Dados do artigo Zhang
lambda <- c(6.04, 8.90)
Q      <- matrix(c(-lambda[1], lambda[1],
                    lambda[2], -lambda[2]), nrow=2,byrow=TRUE)
sigma  <- c(0.44, 0.63) 
r      <- c(1.5, -1.61)
mu     <- r + sigma^2/2
X0     <- 70.5
p0     <- c(0.5,0.5)

h     <- 0.01
beta  <- 10

grid <- expand.grid(seq(-1,1,h),1:2) # malha de estados
G    <- nrow(grid)                   # número de nós na malha
L    <- G/2                          # número de pontos x

g     <- function(x,i)   1-exp(x)
kf    <- function(x,i)   0
b     <- function(x,i)   mu[i]
sig   <- function(x,i)   sigma[i]
q     <- function(x,i,j) Q[i,j]

# Esta função calcula o custo esperado após um
#  passo, seguindo o método da cadeia de Markov
#  aproximada.
expect <- function(n,V)
{
        x <- grid[n,1]
        i <- grid[n,2]
        
        M  <- abs(b(x,i))*h+sig(x,i)^2-q(x,i,i)*h^2
        Dt <- h^2/M
        p_plus  <- (h*max(b(x,i),0) + sig(x,i)^2/2)/M
        p_minus <- (-h*min(b(x,i),0) + sig(x,i)^2/2)/M
        j <- ifelse(i==1,2,1)
        p_j <- q(x,i,j)*Dt
        
        # isso é usado para encontrar a
        #  posição dos vizinhos de x na malha
        bi <- floor((n-1)/L)
        di <- n-bi*L
        
        # calculando a esperança
        e1 <- p_plus *ifelse(di==L,g(x+h),V[di+1+bi*L])
        e2 <- p_minus*ifelse(di==1,g(x-h),V[di-1+bi*L])
        e3 <- p_j*V[n+(-1)^(i-1)*L]
       
        return(exp(-beta*Dt)*(e1 + e2 + e3 + kf(x,i)*Dt))
}

findOptV <- function()
{
    V <- g(grid[,1])
    maxDiff <- 1
    while( maxDiff > 1e-3 )
    {
        maxDiff <- 0
        for( n in 1:G ) # iterando sobre a malha
        {
            x <- grid[n,1]
            i <- grid[n,2]
            Vnew <- min(g(x,i),expect(n,V))
            if( abs((Vnew-V[n])/ifelse(abs(V[n])<0.1,0.1,V[n])) > maxDiff )
                maxDiff = abs((Vnew-V[n])/ifelse(abs(V[n])<0.1,0.1,V[n]))
        
            # atualizando a função valor ótima
            V[n] <- Vnew
        }
        cat(sprintf("Max Diff: %e\n",maxDiff))
     }
    return(V)
}

Vopt <- findOptV()

plot(X0*exp(grid[1:L,1]),g(grid[1:L,1]),lty=2,type="l",ylab="Função Valor",xlab="Preço")
lines(X0*exp(grid[1:L,1]),Vopt[1:L],col="blue",lwd=2)
lines(X0*exp(grid[1:L,1]),Vopt[(L+1):(2*L)],col="red",lwd=2)
lines(X0*exp(grid[1:L,1]),Vopt[1:L]*p0[1]+Vopt[(L+1):(2*L)]*p0[2],col="green",lwd=2)

# Encontrando o momento de parar
S <- rep(0,L)
for( k in 1:L )
{
    x  <- grid[k,1]
    if( g(x,1)*p0[1] + g(x,2)*p0[2] <= expect(k,Vopt)*p0[1] + expect(k+L,Vopt)*p0[2])
         S[k] <- 1
    else
         S[k] <- 0
}
n0 <- min(which(S==0))
n1 <- ifelse(n0==0,n0,n0-1)
n2 <- min(which(S[n0:L]==1))+n1
cat(sprintf("Momento de venda: %f (stop)\n",X0*exp(grid[n1+1,1]-h)))
cat(sprintf("Momento de venda: %f (ganho)\n",X0*exp(grid[n2,1])))
