# Code for Bell (2019). A Measure of Social Coordination and Group Signaling in the Wild

# Define working directory with all the necessary functions or folders
# It should have the following files or folders
# tongatriad2018.rdata; utahtriad2018.rdata; baryplot (folder); modelsim_plotdata.rdata; 
# CodeBell2019C.r; simTestN1000_100c.rdata; restrifit2.rdata; www (folder)
wkdir <- "/Users/pangai/Documents/Dropbox/ActiveProjects/Comparative Classification/PNAS"
# ================================================
# EVOLUTION OF OBJECT ASSOCIATIONS VIA SOCIAL COORDINATION
# ================================================
# define recursion function
dx <- function(x12, x23, d1, d2, d3){

	x13 <- 1 - x12 - x23

	w12 <- x12^2 * (d1 + d2) + x23^2 * (d2) + x13^2 * (d1) + 2*x12*x23 * (d2) + 2*x12*x13 * (d1) + 2*x23*x13 * (0)
	w23 <- x12^2 * (d2) + x23^2 * (d2 + d3) + x13^2 * (d3) + 2*x12*x23 * (d2) + 2*x12*x13 * (0) + 2*x23*x13 * (d3)
	w13 <- x12^2 * (d1) + x23^2 * (d3) + x13^2 * (d1 + d3) + 2*x12*x23 * (0) + 2*x12*x13 * (d1) + 2*x23*x13 * (d3)
			
	wav <- x12 * w12 + x23 * w23 + x13 * w13

	x12p <- x12 * w12 / wav
	x23p <- x23 * w23 / wav
	x13p <- x13 * w13 / wav	
	ans <- list( x12p, x23p, x13p )
	names(ans) <- c("x12p", "x23p", "x13p")
	ans
	 }

dx( x12=0.3, x23=0.3, d1=2, d2=3, d3=1 )

# simulate
covsim <- function(time, x12i, x23i, d1=2, d2=3, d3=1, plotit=TRUE ){

	x12t <- rep(NA, time)
	x23t <- rep(NA, time)
	x13t <- rep(NA, time)

	# starting frequencies
	x12t[1] <- x12i
	x23t[1] <- x23i
	x13t[1] <- 1 - x12i - x23i


	for ( t in 1:(time-1) ){
		dxp <- dx( x12t[t], x23t[t], d1=d1, d2=d2, d3=d3 )
		x12t[t+1] <- dxp$x12p
		x23t[t+1] <- dxp$x23p
		x13t[t+1] <- dxp$x13p
		}
	p <- list( x12t, x23t, x13t )
	names(p) <- c("x12", "x23","x13")	
	
	if( plotit==TRUE ){
		par( mar=c(2,3,0.5,0))
		plot(1:time, x12t, type="l", axes=FALSE, ylim=c(0,1), ylab="", xlab="", col="blue")
		lines( 1:time, x23t, lty=2, col="green" )
		lines( 1:time, x13t, lty=3, col="red" )
		legend( "topright", legend=c("x12", "x23", "x13"), lty=c(1,2,3), box.lwd=0, col=c("blue", "green", "red") )
		axis(side=1)
		axis(side=2)				
		}
	p
}

# x01 x10 common, when will x11 be common?
y <- covsim(time=1000, x12i=0.33333, x23i=0.33333, d1=5, d2=1, d3=1, plotit=TRUE)$x11[1000]
y <- covsim(time=1000, x12i=0.33333, x23i=0.33333, d1=2, d2=2, d3=0.3, plotit=TRUE)$x11[1000]

# understand features of the model
# go over all values of a and d
d1.seq <- seq(0,3,0.05)
d2.seq <- seq(0,3,0.05)
x12eq <- matrix( rep(NA, length(d1.seq)*length(d2.seq) ), nrow=length(d1.seq) )
for ( i in 1:length(d1.seq) ){
	for ( j in 1:length(d2.seq) ) {
		dum <- covsim(time=1000, x12i=0.33333, x23i=0.33333, d1=d1.seq[i], d2=d2.seq[j], d3=0.5, plotit=FALSE)
		x12eq[i,j] <- dum$x12[1000]
	}
}

row.names( x12eq ) <- d1.seq
colnames( x12eq ) <- d2.seq

# ----------
# relationship between coordination and association
par( mar=c(4,4.1,1,1), las=0 )
plot( c(0,2), c(0,2), axes=FALSE, ylab="", xlab="", type="n", lty=2 )
axis(1)#, at=a.seq, labels=a.seq) 
axis(2, las=1)
#lines( c(0,1), c(0,1), lty=2, col="black")
image(x=d1.seq, y=d2.seq, x12eq, col=gray(seq(1,0.5,-0.01), alpha=0.7), add=TRUE )
mtext( expression(paste("Benefits to Coordination on Object 1 (", delta[1], ")",sep="")), side=2, line=3)
mtext( expression(paste("Relative Benefits to Coordination on Object 2 (",delta[2],")")), side=1, line=3)

# ternary plot to understand dynamics
# library(baryplot)
# load baryplot functions
barypath = paste( wkdir,"/baryplot/R/",dir( paste( getwd(),"/baryplot/R", sep="") ), sep="")
lapply( barypath, function(z) source( z ) )

# set payoffs
d1 = 1; d2 = 0.5; d3 = 0.5
	#--------------------------
	# game fitness function
		game = function( p, q, w0=1 ){				
			# fitness for strategies
			W12 = w0 + p * ( d1 + d2 ) + q * d2 + (1 - p - q) * d1
			W23 = w0 + p * d2 + q * (d2 + d3) + ( 1 - p - q) * d3
			W13 = w0 + p * d1 + q * d3 + (1 - p - q) * ( d1 + d3 )
			c( W12, W23, W13 ) # return fitnesses
	
			}		
		#--------------------------
	# ---------------
	# Figure 1
	# ---------------
	# run baryplot
	bary.init()	
	#bary.labels("12","23","13")
	  ## labels for strategies
	  text( 0, 0, "1,3", xpd=NA, adj=c(1,1.5), cex=1.5 ); #left
	  text( 1, 0, "1,2", xpd=NA, adj=c(0,1.5), cex=1.5 ); #right
	  text( 0.5, 1, "2,3", xpd=NA, adj=c(0.5,3), cex=1.5 ); #top
	# phase plot
	bary.phase(thegame = game, length = 0.011)

	# contour plot
	#bary.contour(strat=1, thegame = game)
	#bary.click( thegame=game, arrows=T )
	
#	bary.plotsim( 0.01, 0.99, arrow=TRUE, thegame=game )
	dev <- 0.05
	x12mid <- 1-(2*d1*d2)/(d2*d3+d1*(d2+d3))
	bary.plotsim( x12mid+dev, x23mid-dev/2, arrow=TRUE, thegame=game )
	bary.plotsim( x12mid-dev/2, x23mid+dev, arrow=TRUE, thegame=game )
	bary.plotsim( x12mid-dev/2, x23mid-dev/2, arrow=TRUE, thegame=game )		
	# between 12 and 23
	bary.plotsim(d1/(d1+d3) + dev, d3/(d1+d3) - dev, arrow=TRUE, thegame=game )	
	bary.plotsim(d1/(d1+d3) - dev, d3/(d1+d3) + dev, arrow=TRUE, thegame=game )		# between 23 and 13
	bary.plotsim(0, d2/(d1+d2) - dev, arrow=TRUE, thegame=game )			
	bary.plotsim(0, d2/(d1+d2) + dev, arrow=TRUE, thegame=game )				
	# between 12 and 13
	bary.plotsim(d2/(d2+d3) - dev, 0, arrow=TRUE, thegame=game )			
	bary.plotsim(d2/(d2+d3) + dev, 0, arrow=TRUE, thegame=game )	
	
	bary.point( c(x12mid, x23mid), pch=19)	
	bary.line( c(x12mid, x23mid), c(0.39,0.5+0.45), l=2 )
	eqcex <- 1
	bary.text( c(0.51,0.5+0.52), labels = expression( hat(x)[12]==frac(delta[3](delta[1]+delta[2])-delta[1]*delta[2], delta[2]*delta[3]+delta[1]*(delta[2]+delta[3] ) ) ), lcex=eqcex )			
	bary.text( c(0.57,0.5+0.39), labels = expression( hat(x)[23]==frac(delta[1](delta[2]+delta[3])-delta[2]*delta[3], delta[2]*delta[3]+delta[1]*(delta[2]+delta[3] ) ) ), lcex=eqcex )			
	
	bary.point( c(d3/(d1+d3),d1/(d1+d3)), pch=19)	
	bary.text( c(d3/(d1+d3)+0.12,d1/(d1+d3)), labels = expression( hat(x)[12]==frac(delta[3],delta[1]+delta[3]) ), lcex=eqcex )					

	bary.point( c(0,d1/(d1+d2)), pch=19)	
	bary.text( c(0-0.12,d1/(d1+d2)), labels = expression( hat(x)[23]==frac(delta[1],delta[1]+delta[2]) ), lcex=eqcex )			
	
	bary.point( c(d3/(d2+d3),0), pch=19)					
	bary.text( c(d3/(d2+d3),0-0.1), labels = expression( hat(x)[12]==frac(delta[3],delta[2]+delta[3]) ), lcex=eqcex )	
	
# ================================================
# STATISTICAL MODEL 
# ================================================
# Model with function as a competition between two forces for classification:
# Social Coordination and Inherent Similarity

# ---------------------------------------------
# GET TONGA DATA 
# ---------------------------------------------
# Loads the classification data from Tonga
load(paste0(wkdir,"/tongatriad2018.rdata"))

# ---------------------------------------------
# GET UTAH DATA 
# ---------------------------------------------
# Loads the classification data from Utah
load(paste0(wkdir,"/utahtriad2018.rdata"))

# DATA LOADING COMPLETE WITH TWO OBJECTS
# Tongan data: d.tonga
# Utah data: d.not

# -----------------------------
# STATISTICAL MODEL OF SOCIAL COORDINATION
# exploring key features
# -----------------------------
# define inverse logit function
logit <- function(x) exp(x)/(1+exp(x)) 

# define function with recursions
xtl <- function( x12,x13,d1,d2,d3,Bs=1 ){
	x23 <- 1 - x12 - x13

	x12p <- 2*x12*x23*(1/2)*logit(Bs*d2) + 2*x12*x13*(1/2)*logit(Bs*d1) + 2*x23*x13*(1-logit(Bs*d3)) + x12^2*(1 - (1/2)*logit(Bs*(d1-d2)) - (1/2)*logit(Bs*(d2-d1)) ) + x23^2*(1/2)*logit(Bs*(d2-d3)) + x13^2*(1/2)*logit(Bs*(d1-d3))  
	
	x13p <-  2*x12*x23*(1 - logit(Bs*d2)) + 2*x12*x13*(1/2)*logit(Bs*d1) + 2*x23*x13*(1/2)*logit(Bs*d3) + x12^2*(1/2)*logit(Bs*(d1-d2)) + x23^2*(1/2)*logit(Bs*(d3-d2)) + x13^2*(1 - (1/2)*logit(Bs*(d1-d3))-(1/2)*logit(Bs*(d3-d1))) 
	
	x23p <- 1 - x12p - x13p
	ans <- list( x12p, x13p, x23p )
	names(ans) <- c("x12", "x13", "x23")
	ans
}

xtl(x12=0.3,x13=0.3,d1=1,d2=1,d3=0.3,Bs=0.1)

# function to run and/or plot evolutionary
# dynamics of association
assoc.sim.l <- function(time=100, x12i=0.3,x13i=0.3,d1=0.2,d2=0.3,d3=0.4, Bs=1/3, plotit=FALSE ){

	x12.tm <- rep(NA, time)
	x13.tm <- rep(NA, time)
	x23.tm <- rep(NA, time)

	# starting frequencies
	x12.tm[1] <- x12i
	x13.tm[1] <- x13i
	x23.tm[1] <- 1-x12i-x13i	

	for ( t in 1:(time-1) ){
		xtp <- xtl( x12=x12.tm[t], x13=x13.tm[t],d1=d1,d2=d2,d3=d3,Bs=Bs )
		x12.tm[t+1] <- xtp$x12
		x13.tm[t+1] <- xtp$x13
		x23.tm[t+1] <- xtp$x23		
		}
	x <- list( x12.tm, x13.tm, x23.tm )
	names(x) <- c("x12", "x13","x23")	

	if( plotit==TRUE ){
		par( mar=c(2,3,0.5,0))
		plot(1:time, x12.tm, type="l", axes=FALSE, ylim=c(0,1), ylab="", xlab="", col="blue")
		lines( 1:time, x13.tm, lty=2, col="red" )
		lines( 1:time, x23.tm, lty=3, col="green" )
		legend( "topright", legend=c("x12", "x13", "x23"), lty=c(1,2,3), box.lwd=0, col=c("blue", "red", "green") )
		axis(side=1)
		axis(side=2)		
		}
	x
}

# check model simulation
y <- assoc.sim.l(time=100, x12i=0.3,x13i=0.3,d1=2,d2=8, d3=10, Bs=0.3, plotit=TRUE)
yt <- data.frame(y$x12,y$x13,y$x23)
rowSums(yt)
yt[100,]

# use endpoint of simulation as
# a theoretical prediction to 
# construct Pr( x | d1, d2, d3, Bs )

# -----------------
# Describe nature of the stats model with 
# a few key figures
# -----------------
# Note LONG TIME! note these for loops take about a 1/2 hour 
# alternatively load this data and go straight to plotting below
load(paste0(wkdir,"/modelsim_plotdata.rdata"))

d2seq <- seq(0,10,0.01)
d1seq <- seq(0,10,0.01)
y1 <- matrix( rep(NA,length(d2seq)*length(d1seq)), nrow=length(d2seq) )
for( i in 1:length(d2seq) ){
	for( j in 1:length(d1seq) ){
	y1[i,j] <- assoc.sim.l(time=100, x12i=0.3, x13i=0.3, d1=d1seq[j], d2=d2seq[i], d3=0.3, Bs=0.3, plotit=FALSE)$x12[100]
			}
	}	

# This is Figure 3
pdf(file="modelsim_plot.pdf")
quartz(width=11, height=5)
par( mar=c(4,4.2,1,1), mfrow=c(1,1) )
plot( c(0,max(d2seq)), c(0,max(d1seq)), axes=FALSE, ylab="", xlab="", type="n", lty=2 )
axis(1) 
axis(2)
image(x=d2seq, y=d1seq, z=y1, col=gray(seq(1,0,-0.01), alpha=1), add=TRUE )
mtext( expression(paste("Benefits to Coordination on Item 1 ", (delta[1]))), side=2, line=3, cex=0.8)
mtext( expression(paste("Benefits to Coordination on Item 2 ", (delta[2]))), side=1, line=3, cex=0.8)
contour(x=d2seq, y=d1seq, y1, add=TRUE )
dev.off()

# -----------------
# do initial frequencies matter?
# the simulation below suggets no
freqseq <- runif(100,0, 0.7)
simInits <- sapply( freqseq, function(z) {
	y <- assoc.sim.l(time=100, x12i=z,x13i=0.3,d1=0.7,d2=0.3,d3=0.4, Bs=1/3, plotit=FALSE)
	yt <- data.frame(y$x12,y$x13,y$x23)
	#rowSums(yt)
	yt[100,]
	})
plot( x=c(0,1), y=c(0,1), type="n", xlab="Initial freq. of x12", ylab="end freqs. of simulation")
for( i in 1:dim(simInits)[2] ) points( x=rep(freqseq[i],3), y=simInits[,i], lty=1:3 )

# how do equilibria respond to the parameter space
# vary one parameter while holding the rest constant

sim.par <- function( wpar, parseq ){
	ans <- sapply( parseq, function(z){
		pargs <- list( z )
		names(pargs) <- wpar		
		#print(pargs)
		y <- do.call("assoc.sim.l", args=pargs)
		c(y$x12[100],y$x13[100],y$x23[100])
		})	
	}

plot.parsim <- function(x, parseq){
		par( mar=c(2,3,0.5,0))
		plot(parseq, x[1,], axes=FALSE, ylim=c(0,1), type="n", ylab="", xlab="" )
		lines( parseq, x[1,], lty=1, col="blue" )
		lines( parseq, x[2,], lty=2, col="red" )
		lines( parseq, x[3,], lty=3, col="green" )
		legend( "topright", legend=c("x12", "x13", "x23"), lty=c(1,2,3), box.lwd=0, col=c("blue", "red", "green") )
		axis(side=1)
		axis(side=2)		
	}

# d1, d2, d3, d4, d5, d6
# selection weight Bs, and u

sim.Bs <- sim.par( wpar="Bs", parseq=seq(0,10,0.01))
plot.parsim( x=sim.Bs, parseq=seq(0,10,0.01))

sim.d1 <- sim.par( wpar="d1", parseq=seq(0,10,0.01))
plot.parsim( x=sim.d1, parseq=seq(0,10,0.01))

sim.d2 <- sim.par( wpar="d2", parseq=seq(0,10,0.01))
plot.parsim( x=sim.d2, parseq=seq(0,10,0.01))

sim.d3 <- sim.par( wpar="d3", parseq=seq(0,10,0.01))
plot.parsim( x=sim.d3, parseq=seq(0,10,0.01))

# -----------------------------
# Model of inherent similarity
# -----------------------------
# Train model with data from non-Tongans

# function to count responses to each triad
triadCount <- function( y=y, triad=triad, combin )
		{
	# needs data triad (number combination)
	# and y, which is a matrix of zeros and ones
	# x are the parameters to the model
	ntri <- dim(triad)[1]
	combin <- t(combn(1:6,3))
	count <- matrix( rep(0, ntri*3), nrow=ntri )

	pr12 <- pr13 <- pr23 <- rep(NA, ntri)	
	
	for ( i in 1:ntri ){
		pr <- triad[i,]

		# which row of parameters to use
		wrow <- which( sapply( 1:ntri, function(z) any(combin[z,1:3]==pr[1]) & any(combin[z,1:3]==pr[2]) & any(combin[z,1:3]==pr[3] ) )	)	
		
		# for each triad, count
		add.index <- which( combin[wrow,] == triad[i, y[i,2:4]==1 ] )
		count[wrow,add.index] <- count[wrow,add.index] + 1		
		}
	count							
	}

# test
triadtest <- d.not[[1]]$motif_tricomb[[1]]
ytest <- d.not[[1]]$motif_triad_data[[1]]

triadCount( y = ytest, triad = triadtest, combin = combin )


inherent.similarity <- function(d){
	combin <- t(combn(1:6,3))
	ans <- matrix( rep(0,dim(combin)[1]*dim(combin)[2] ), nrow=dim(combin)[1] )
	for ( i in 1:length(d) ){
		triad <- d[[i]]$motif_tricomb[[1]]
		y <- d[[i]]$motif_triad_data[[1]]
		ans <- ans + triadCount( y = y, triad = triad, combin = combin ) 
		}
	ans	
}

inher.sim.freq <- inherent.similarity( d.not ) / length( d.not )

################################
# Likelihood model

mtriadProb <- function( y=y, triad=triad, delta, Bs.par, u.par, inher.sim )
		{
	# first check conditional probabilities
	# make sense, if not then return low likelihood
						
	# needs data triad (number combination)
	# and y, which is a matrix of zeros and ones
	# x are the parameters to the model
	combin.motif <- t(combn(1:6,3))
	ntri <- dim(combin.motif)[1] # no. of triads
	
	pr12 <- pr13 <- pr23 <- rep(NA, ntri)	
	
	llik <- 0 # set up likelihood	
	for ( i in 1:ntri ){
		pr <- triad[i,]
		# what should the starting frequencies be of the
		# evolutionary model? from simulations above, the 
		# equilibrium is not sensitive to initial frequencies
		# specify near-equal frequencies of the three strategies
		em <- 	assoc.sim.l(time=100, x12=0.333,x13=0.333,d1=delta[pr[1]],d2=delta[pr[2]],d3=delta[pr[3]],Bs=Bs.par, plotit=FALSE)
		emt <- data.frame(em$x12,em$x13,em$x23)
		# check results
		# rowSums(emt)

		# -----------
		# trained model for inherent similarity
		isr <- which( sapply( 1:ntri, function(z) any(combin.motif[z,1:3]==pr[1]) & any(combin.motif[z,1:3]==pr[2]) & any(combin.motif[z,1:3]==pr[3] ) )	)	

		# prob of pair 1,2
		pr1 <- u.par * emt[100,1] + ( 1 - u.par ) * inher.sim[isr,3]
		# prob of pair 1,3
		pr2 <- u.par * emt[100,2] + ( 1 - u.par ) * inher.sim[isr,2]
		# prob of pair 2,3
		pr3 <- u.par * emt[100,3] + ( 1 - u.par ) * inher.sim[isr,1]

		# final log probabilities
		# first and second object grouped
		pr12[i] <- log( pr1 * ( 1 - pr2 ) * ( 1 - pr3 ) )
		# first and third object grouped		
		pr13[i] <- log( pr2 * ( 1 - pr1 ) * ( 1 - pr3 ) )
		# second and third object grouped		
		pr23[i] <- log( pr3 * ( 1 - pr1 ) * ( 1 - pr2 ) )
		
		obj <- which( y[i,2:4] == 1 )
		# probability
		prob <- ifelse( obj==1, pr23[i], ifelse( obj==2, pr13[i], ifelse( obj==3, pr12[i], 0 ) ) )
		llik <- llik + prob # add the log probabilities
		}
	-llik	
	}

m1loglik <- function( x ){
	# all values are >zero with
	# u from zero to one
	# so use logit transformation
	
	nObj <- 6 # six objects
	int <- t(combn(1:nObj,3)) # triad combinations
	ntri <- dim(int)[1] # no. of triads

	# second set of parameters in x, delta,
	# benefits to social coordination on one
	# signal is passed to the evolutionary model
	delta <- exp(x[1:nObj]) 
	
	# parameter scaling payoff differentials
	# to strength of selection
	Bs.par <- exp(x[(nObj+1)])
	
	# parameter weighing the relative weight of
	# social coordination explaining data versus
	# inherent similarity
	u.par <- plogis(x[(nObj+2)])

	# number of participants
	# Tongan data: d.tonga
	nObs <- length(d.tonga)

	neglogLik <- rep(0,nObs)
	for( i in 1:nObs ){
		neglogLik[i] <- mtriadProb( y=d.tonga[[i]]$triad_data, triad=d.tonga[[i]]$tricomb, delta=delta, Bs.par=Bs.par, u.par=u.par, inher.sim=inher.sim.freq  )
		}
	sum(neglogLik)		
}

# test likelihood function
# initial parameter values into vector  x
delta.inits <- runif( 6)
Bs.init <- runif( 1 )
u.init <- runif( 1 )
x.inits <- c( delta.inits, Bs.init, u.init )
m1loglik (x.inits)
# very slow calculation

# load likelihood coded in C 
source(paste0(wkdir,"/CodeBell2019C.r"))

################################
# INVESTIGATE THE SURACE OF THE LIKELIHOOD
# FOR TWO PARAMETERS

# prep data for C function of negLogLikelihood
y <- d.tonga[[1]]$triad_data[,2:4]
for( i in 2:length(d.tonga) ) y <- cbind( y,d.tonga[[i]]$triad_data[,2:4] )
y <- as.integer(y)
alphacomb <- as.integer(t(combn(1:6,2))) # alpha_ij from 15x2 matrix to vector
combin <- t(combn(1:6,3))
tricomb <- d.tonga[[1]]$tricomb
for( i in 2:length(d.tonga) ) tricomb <- cbind( tricomb,d.tonga[[i]]$tricomb )
tricomb <- as.integer(tricomb)

# loop over values of u and Bs
zseq <- seq( 0, 1, .1 ) # range for u
yseq <- seq( 0, 1, .1 ) # range for Bs
mB.surf <- matrix( rep(NA,length(zseq)*length(yseq)), nrow=length(zseq) )
for( i in 1:length(zseq) ){
	for(j in 1:length(yseq) ){
		x.pars <- c( exp(x.inits[1:6]), plogis( c(zseq[i],yseq[j] )  ) )
		mB.surf[i,j] <- -funLL( x = x.pars, y=y, alphacomb=alphacomb, nObj=as.integer(6), combin=as.integer(combin), tricomb=tricomb, inhersim=as.numeric(inher.sim.freq), candLL=numeric(1), N=as.integer(length(d.tonga)), nTri=as.integer(20), debug=numeric(1) )$candLL
		}		
	}

library(graphics)
# find point where the maximum likelihood is
minval <- which( mB.surf==min(mB.surf,na.rm=TRUE), arr.ind=TRUE )
zseq[minval[1]]
yseq[minval[2]]
# plot the contour with a red point showing the MLEs
quartz()
mB.surf2 <- ifelse( mB.surf>10000, NA, mB.surf )
lvls <- 13
filled.contour(x=zseq, y=yseq, z=mB.surf2, xlab="u", ylab="Bs", nlevels=lvls, col=gray(0:lvls/lvls), plot.axis={ 	points(zseq[minval[1]],yseq[minval[2]],pch=3, col="red", bg="black", lwd=2)	} )

# likelihood looks smooth
# no discontinuities, yeah!

################################
# DATA SIMULATION AND TEST OF 
# STATISTICAL ROUTINE
################################
# delta: signaling value
# d1, d2, d3, d4, d5, d6
# selection weight Bs and weight u

# add two more parameters for a total of 8

# data simulation function
data.sim <- function( x=x, inher.sim, N=100 ){
	
	nObj <- 6 # six objects
	# needs data triad (number combination)
	combin.motif <- t(combn(1:nObj,3))
	ntri <- dim(combin.motif)[1] # no. of triads

	# first set of parameters in x: delta
	# benefits to social coordination on one
	# signal is passed to the evolutionary model
	delta <- x[1:nObj]
	
	# parameter scaling payoff differentials
	# to selection
	Bs.par <- x[(nObj+1)]	# with constraint
	
	# parameter weighing the relative weight of
	# social coordination explaining data versus
	# inherent similarity
	u.par <- x[(nObj+2)]			
	
	# use lapply to simulate observations
	# across N individuals
	dsim <- lapply( 1:N, function(z){
		sim.triad <- matrix( rep(NA,ntri*4),nrow=ntri)	
		for ( i in 1:ntri ){
		pr <- combin.motif[i,]
		# calculate probability of each pairing
		# given a triad
		
		# for each triad, three models for three outcomes
		
		# what should the starting frequencies be of the
		# evolutionary model? from simulations above, the 
		# equilibrium is not sensitive to initial frequencies
		# specify near-equal frequencies of the three strategies
		em <- 	assoc.sim.l(time=100, x12=0.333,x13=0.333,d1=delta[pr[1]],d2=delta[pr[2]],d3=delta[pr[3]],Bs=Bs.par, plotit=FALSE)
		emt <- data.frame(em$x12,em$x13,em$x23)

		# -----------
		# trained model for inherent similarity
		isr <- which( sapply( 1:ntri, function(z) any(combin.motif[z,1:3]==pr[1]) & any(combin.motif[z,1:3]==pr[2]) & any(combin.motif[z,1:3]==pr[3] ) )	)	

		# prob of pair 1,2
		pr1 <- u.par * emt[100,1] + ( 1 - u.par ) * inher.sim[isr,3]
		# prob of pair 1,3
		pr2 <- u.par * emt[100,2] + ( 1 - u.par ) * inher.sim[isr,2]
		# prob of pair 2,3
		pr3 <- u.par * emt[100,3] + ( 1 - u.par ) * inher.sim[isr,1]

		# with probabilites derived from the model
		# draw random number to select outcome
		rn <- runif( 1 )
		if( rn<=pr1 ) sim.triad[i,] <-c(i,0,0,1)
		if( rn>pr1 & rn<=(pr1+pr2) ) sim.triad[i,] <-c(i,0,1,0)
		if( rn>(pr1+pr2) ) sim.triad[i,] <- c(i,1,0,0) 
			} # close for loop
			sim.triad
		}) # end of lapply
		dsim
	}

# run estimation routine
nRuns <- 1000 # number of simulations
nObj <- 6 # number of objects
pars.est <- pars.sim <- matrix( rep(NA,nRuns*24), nrow=nRuns )
llikdiff <- rep(NA, nRuns)
# set up the combination
# unlike the data, it is the same
# order in the simulation
tricomb.sim <- t(combn(1:nObj,3))
for( i in 2:100 ) tricomb.sim <- cbind( tricomb.sim, combin )
tricomb.sim <- as.integer(tricomb.sim)

# THIS TAKES A WHILE (1 day at least)
# YOU CAN LOAD THIS DATA INSTEAD
# AND SKIP TO THE PLOTS BELOW
load("simTestN1000_100c.rdata")

for ( i in 1:nRuns ){
delta.inits <- runif( 6, 0, 6 )
Bs.init <- runif( 1, 0.005, 3 )
u.init <- runif( 1, 0, 1 )
x.inits <- c( delta.inits, Bs.init, u.init )	
#x.inits.ql <- qlogis( x.inits ) # logit function since inverse logit used in routine
pars.sim[i,] <- x.inits
d.sim <- data.sim( x.inits, inher.sim =inher.sim.freq, N=100 )

y.sim <- d.sim[[1]][,2:4]
for( j in 2:length(d.sim) ) y.sim <- cbind( y.sim,d.sim[[j]][,2:4] )
y.sim <- as.integer(y.sim)

mfit.tri.wrp.sim <- function(xinits) {
	inis <-  as.numeric( c( exp(xinits[1:7]), plogis(xinits[8]) ) )
	-funLL( x = inis, y=y.sim, alphacomb=alphacomb, nObj=as.integer(6), 	combin=as.integer(combin), tricomb=tricomb.sim, inhersim=as.numeric(inher.sim.freq), candLL=numeric(1), N=as.integer(length(d.sim)), nTri=as.integer(20), debug=numeric(1) )$candLL
}

pert.xinits <- c(log(abs(x.inits[1:7]+runif(7,-0.1,0.1))),qlogis(abs(min(c(x.inits[8]+runif(1,-0.1,0.1),1) ) ) ) )

check <- mfit.tri.wrp.sim (xinits=pert.xinits )

if( !is.na(check) ){

	tryCatch({
    	print(i)
		ans <- optim( pert.xinits, mfit.tri.wrp.sim, method="Nelder-Mead", hessian=FALSE, control=list(trace=TRUE) )
		llikdiff[i] <- mfit.tri.wrp.sim(c(log(x.inits[1:7]),qlogis(x.inits[8]) ) ) - ans$val 
		pars.est[i,]<-c( exp(ans$par[1:7]), plogis(ans$par[8]) ) 	    
	  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

	}
}

# parameter labels
plabs <- c( sapply(1:6, function(z) bquote( delta[.(z)] ) ),expression(beta[S]), expression(italic(u)) )

# supplementary figures to
# understand the performance of 
# the MLE approach
png( file="simParCompareReduced.png" , width=700, height=500, res=100)
par(mar=c(2,2,1,1), mfrow=c(4,4), cex.lab=0.8, las=1)
lms <- c(7,7,7,7,7,7,4,1)
for( j in 1:8 ){
lims <- min( 8, max( c( pars.sim[1:nRuns,j],pars.est[1:nRuns,j]), na.rm=T ) )
plot( x=c(0,lms[j]), y=c(0,lms[j]), xlim=c(0, lms[j] ), ylim=c(0,lms[j]), type="l", lty=2, ylab="", xlab="", axes=FALSE)
axis(1)
axis(2)
points( x=pars.sim[1:nRuns,j], y=pars.est[1:nRuns,j], col="darkgray", cex=0.8, xpd=TRUE )
text(x=(1/10)*lms[j],y=(5/6)*lms[j], labels=plabs[j], cex=1.2, col="black")
	}

for( j in 1:8 ){
	pardiff <- pars.est[1:nRuns,j]-pars.sim[1:nRuns,j]
hist( pardiff[-5<=pardiff & pardiff <=5 & !is.na(pardiff)], main="", col="gray", xlim=c(-lms[j],lms[j]) )
av <- format( mean( pardiff, na.rm=TRUE  ), digits=2)
med <- format( median( pardiff, na.rm=TRUE  ), digits=2)
legend( "topright", legend=plabs[j], box.lwd=0 , bg=gray(1, alpha=0))
legend( "right", legend=c(paste0(" mean = ",av),paste0(" median = ",med)), box.lwd=0 , bg=gray(1, alpha=0), yjust=2, cex=0.8)
	}
dev.off()

################################
# now that the ML estimates 
# appear unbiased, use
# an optimization algorithm
# to find the estimates at the 
# minimum negative log-likelihood

mfit.tri.wrp <- function(xinits) {
	inits <-  as.numeric( c( exp(xinits[1:7]), plogis(xinits[8]) ) )
	-funLL( x = inits, y=y, alphacomb=alphacomb, nObj=as.integer(6), 	combin=as.integer(combin), tricomb=tricomb, inhersim=as.numeric(inher.sim.freq), candLL=numeric(1), N=as.integer(length(d.tonga)), nTri=as.integer(20), debug=numeric(1) )$candLL
}

mfit.tri.wrp(x.pars)
mfit <- optim( x.pars, mfit.tri.wrp, method="Nelder-Mead", hessian=TRUE, control=list(trace=TRUE, maxit=100))
exp(mfit$par[1:7])
plogis(mfit$par[8])
sqrt( diag( solve( mfit$hessian ) ) )

# point estimate precision is 
# hard to estimate 

# try Monte Carlo methods
# set up the input for the 
# Metropolis-Hastings sampler
# use functions loaded from 
# CodeBell2019C

# parameter labels
par.labs <- c(paste0("d",1:6), sapply( 1:15, function(z) paste0( "a", paste0(c(combn(1:6,2)[,z]), collapse="" ) ) ), "Bs", "Bd", "u" )

# classification data
y <- d.tonga[[1]]$triad_data[,2:4]
for( i in 2:length(d.tonga) ) y <- cbind( y,d.tonga[[i]]$triad_data[,2:4] )
y <- as.integer(y)

# parameter alpha subscripts
alphacomb <- as.integer(t(combn(1:6,2)))

# combinations of 3 out of 6
combin <- t(combn(1:6,3))

# combinations data
tricomb <- d.tonga[[1]]$tricomb
for( i in 2:length(d.tonga) ) tricomb <- cbind( tricomb,d.tonga[[i]]$tricomb )
tricomb <- as.integer(tricomb)

# sampler specifications
nsamples=1000000 # how many samples

# step sizes for each parameter
#eps = c( 0.05, 0.06, 0.07, 0.12, 0.05, 0.05, 0.1, 0.09, 0.06, 0.06, 0.06,     0.1, 0.06, 0.07, 0.07, 0.08, 0.07, 0.1, 0.09, 0.08, 0.08, 0.1, 0, 0.05)
eps = c( 0.05, 0.06, 0.07, 0.3, 0.05, 0.05, rep(0,15), 0, 0, 0.05)
#eps = c( 0.1, 0.1, 0.1, 0.07, 0.1, 0.1, 0.1, 0.1, 0.07, 0.06, 0.07, 0.08, 0.06, 0.07, 0.07, 0.08, 0.12, 0.13, 0.1, 0.2, 0.12, 0, 0.11, 0.06)
#eps <- ifelse( is.na(se.new) | se.new > 10, 1/100 , se.new/100 )

# pre-allocate jumps across with each sample
npar=as.integer(length(par.labs))
jump = runif( npar*nsamples, -rep(eps,nsamples), rep(eps, nsamples) ) 

# starting parameter values
mh.inits <- as.numeric( c(1.3680970267426, 1.43164564094116, 1.6005294831759, 5.30864212886106, 0.833522188662519, 1.1840217160905, mfit$par[7:21], -0.915401, -1000, 2.22356612462538 ) ) #-0.1478382, 2.135469) )#+runif(24,0,0) )

# test likelihood with initial parameters
mfit.tri.wrp( mh.inits )
mfit.tri.wrp( mfit$par )

source(paste0(wkdir,"/CodeBell2019C.r"))

# MCMC SAMPLING TAKES ABOUT 2 DAYS
# alternatively load this data
# and skip to plotting
load(paste0(wkdir,"/restrifit2.rdata"))

mfit.tri = funMH( npar=npar, m=nsamples, jump=jump, unifs=runif(nsamples), out=c(mh.inits, rep(0,npar*(nsamples-1))), x = mh.inits, y=y, alphacomb=alphacomb, nObj=as.integer(6), combin=as.integer(combin), tricomb=tricomb, inhersim=as.numeric(inher.sim.freq),  N=as.integer(length(mres6)), nTri=as.integer(20), debug=numeric(1), aLL=numeric(nsamples) )	
organize( mfit=mfit.tri, model="trifit2",  inits=mfit.tri$x, par.labels=par.labs, sam=nsamples )

# FIGURE 2
transf <- function(x) c(exp(x[1:6]),plogis(x[24]))
est <- transf( colMeans( mfit$samples[500000:1000000,]) )
est.lb <- transf( sapply( 1:24, function(z) quantile( mfit$samples[500000:1000000,z], probs=0.025) ) )
est.ub <- transf( sapply( 1:24, function(z) quantile( mfit$samples[500000:1000000,z], probs=0.975) ) )
data.frame( est, est.lb, est.ub )

imageDir <- paste0(wkdir,"/www/")
pdf( file=paste0(wkdir,"/probChoice.pdf") )
quartz(width=6.5, height=4.7)
par(mfrow=c(1,1))
delt <- est[1:6]
delt.lb <- est.lb[1:6]
delt.ub <- est.ub[1:6]
b <- 0.25 # image spacer
delt.vals <- delt
delt.or <- order(delt.vals)
par(las=1, mar=c(5,6,1,6), cex=0.8)
mv <- 0.08
plot( c(1,6), c(0,7), type="n", axes=FALSE, ylab="",xlab="" )
points( 1:6-mv, delt.vals[delt.or], pch=16, cex=2)
axis(1)
axis(2)

library(jpeg)
for ( i in 1:nObj ){
whr <- delt.or[i]
img<-readJPEG(paste0(imageDir,type="m",whr,".jpeg"))
rasterImage(img, xleft=i-b,ybottom=-1.5,xright=i+b,ytop=-0.5, xpd=TRUE)
lines( c(i-mv,i-mv), c(delt.lb[whr],delt.ub[whr]) )
lines( c(i-0.1-mv,i+0.1-mv), c(delt.lb[whr],delt.lb[whr])  )
lines( c(i-0.1-mv,i+0.1-mv), c(delt.ub[whr],delt.ub[whr])  )
}

mtext(side=2, expression( paste("Signaling value of motif, ",italic(delta[i]), sep="")), line=4.5, las=0)
utext <- bquote(italic(u) == .(format(est[7],digits=2 )) )	
utext <- bquote(italic(u) == .(paste(format(est[7],digits=3 )," (",format(est.lb[7],digits=3 ),", ",format(est.ub[7],digits=3 ),")",sep="") ))	
text( x=1, y=6.5, utext, adj=0, cex=1.2 )

# coefficient of selection Beta_s
Bs.set <- 0.2858959
Bs.set * delt.vals
invlgt.Bs <- plogis( Bs.set * delt.vals )
invlgt.Bs.lb <- plogis( Bs.set * delt.lb )
invlgt.Bs.ub <- plogis( Bs.set * delt.ub )
#plot(logit(Bs.set * delt.vals[delt.or]))

# add right axis
minaxis <- plogis(Bs.set*0)
maxaxis <- plogis(Bs.set*7)
intvls <- ( maxaxis - minaxis ) / 7
lbsl <- seq( minaxis, maxaxis, intvls)
axis(4, at=0:7, labels=format(lbsl,digits=2), las=1)
mtext(side=4, expression( paste("Probability of associations with motif, H","(",italic( beta[s] ),italic(delta[i]),")", sep="")), line=4.5, las=0)
scl <- function( x ) ( x - minaxis ) / intvls # scale points
for ( i in 1:nObj ){
whr <- delt.or[i]
points( i+mv, scl( invlgt.Bs[whr] ), pch=2, cex=1 )
points( i+mv, scl( invlgt.Bs[whr] ), pch=2, cex=2 )

lines( c(i+mv,i+mv), c(scl(invlgt.Bs.lb[whr]),scl(invlgt.Bs.ub[whr]) ) )
lines( c(i+mv-0.1,i+mv+0.1), c(scl(invlgt.Bs.lb[whr]),scl(invlgt.Bs.lb[whr]) )  )
lines( c(i+mv-0.1,i+mv+0.1), c(scl(invlgt.Bs.ub[whr]),scl(invlgt.Bs.ub[whr]) )  )
}

dev.off()