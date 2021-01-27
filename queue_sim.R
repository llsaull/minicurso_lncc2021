# Simulações de caminhos do modelo de fila apresentado
#  para ilustrar aproximações por tráfego pesado.
library(latex2exp)
T = 1
n = 10
ld = 1
la = 0.9*n
rho = la/(n*ld)

X = 0
t = 0
tc= 0
h = 0.001
while( tc <= T){
    Dt = rexp(1,la)
    W  = rexp(1,ld)
    Xn = Xn + W
    hc = 0
    while(hc < Dt){
        hc = hc + h
        tc = tc + h
        Xn = max(Xn-n*h,0)
        X = c(X,Xn)
        t = c(t,tc)
    }
    hc = Dt-hc
    tc = tc + hc
    Xn = max(Xn-n*hc,0)
    X = c(X,Xn)
    t = c(t,tc)
    cat(sprintf("time %f\n",tc))
}

plot(t,X/sqrt(n),type="l",col = "blue", ylab=TeX("$X/\\sqrt{n}$"))

