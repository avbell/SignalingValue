# Metropolis-Hastings and other computations via C is faster!

require(inline) # wrapper for C code
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# SIMULATION OF EVOLUTIONARY MODEL
evolSimL <- "
	int time=100; // number of timesteps
	double d1d = *d1, d2d = *d2, d3d = *d3;
	double Bsd = *Bs;	
	double x12[time], x13[time], x23, x12d, x13d;
	x12[0] = 1.0/3.0;
	x13[0] = 1.0/3.0;
	// loop over time
		for( int t=1; t<time+1; t++ ){
			x12d = x12[t-1];
			x13d = x13[t-1];
			x23 = 1.0 - x12d - x13d;						
			x12[t] = 2.0 * x12d * x23 * ( 1.0 / 2.0 ) * pow( exp(1), Bsd * d2d )/( 1.0 + pow( exp(1), Bsd * d2d ) ) + 2.0 * x12d * x13d * ( 1.0 / 2.0 ) * pow( exp(1), Bsd * d1d )/( 1.0 + pow( exp(1), Bsd * d1d ) ) + 2.0 * x23 * x13d * ( 1.0 - pow( exp(1), Bsd * d3d ) / ( 1.0 + pow( exp(1), Bsd * d3d) ) ) + pow(x12d,2.0) * ( 1.0 - ( 1.0 / 2.0 ) * ( pow( exp(1), Bsd * ( d1d - d2d ) ) / ( 1.0 + pow( exp(1), Bsd * ( d1d - d2d ) ) ) ) - ( 1.0 / 2.0 ) * ( pow( exp(1), Bsd * ( d2d - d1d ) ) / ( 1.0 + pow( exp(1), Bsd * ( d2d - d1d ) ) ) ) ) + pow(x23,2.0) * ( 1.0 / 2.0 ) * ( pow( exp(1), Bsd * ( d2d - d3d ) )/ ( 1.0 + pow( exp(1), Bsd * ( d2d - d3d ) ) ) ) + pow(x13d,2.0) * ( 1.0 / 2.0 ) * ( pow( exp(1), Bsd * ( d1d - d3d ) )/( 1.0 + pow( exp(1), Bsd * ( d1d - d3d ) ) ) );
									
			x13[t] = 2.0 * x12d * x23 * ( 1.0 - ( pow( exp(1), Bsd * d2d )/(1.0 + pow( exp(1), Bsd * d2d ) ) ) ) + 2.0 * x12d * x13d * ( 1.0 / 2.0 ) * ( pow( exp(1), Bsd * d1d ) / ( 1.0 + pow( exp(1), Bsd * d1d ) ) )  + 2.0 * x23 * x13d * ( 1.0/2.0 ) * ( pow( exp(1), Bsd * d3d ) / ( 1.0 + pow( exp(1), Bsd * d3d ) ) ) + pow(x12d,2.0) * ( 1.0 / 2.0 ) * ( pow( exp(1), Bsd * ( d1d - d2d )  )/( 1.0 + pow( exp(1), Bsd * ( d1d - d2d ) ) ) ) + pow(x23,2.0) * ( 1.0 / 2.0 ) * ( pow( exp(1), Bsd * ( d3d - d2d ) )/( 1.0 + pow( exp(1), Bsd * ( d3d - d2d ) ) ) ) + pow(x13d,2.0) * ( 1.0 - ( 1.0 / 2.0 ) * ( pow( exp(1), Bsd * (d1d - d3d ) ) / ( 1.0 + pow( exp(1), Bsd * ( d1d - d3d ) ) ) ) - ( 1.0 / 2.0 ) * ( pow( exp(1), Bsd * (d3d - d1d ) ) / ( 1.0 + pow( exp(1), Bsd * (d3d - d1d ) ) ) ) );			
					}
			*x12eq = x12[time];
			*x13eq = x13[time];
	"
calevolSimL <- signature( x12eq="numeric", x13eq="numeric", d1="numeric", d2="numeric", d3="numeric", Bs="numeric" )


# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# TRIAD CLASSIFICATION MODEL LIKELIHOOD

LLtriad <- "
	// using tricomb, or various combinations in the data
	double Bs = x[*nObj]; // strength of selection on singles
	double upar = x[*nObj+1]; // weight on social influence		
	double pr12,pr13,pr23,lpr12f,lpr13f,lpr23f,x12eq=0,x13eq=0; // parameter values
	double d1,d2,d3;
	int cmb1,cmb2,cmb3;
	int f1r=0,f1c=0,f2r=0,f2c=0,f3r=0,f3c=0;	 // for search function
	// CALCULATE LOG-LIKELIHOOD
		// loop over observations (individuals) and triads
		for( int i=0; i < *N; i++ ){
			for( int j=0; j < *nTri; j++){
				cmb1 = tricomb[i*60 + 0 * *nTri + j]; // user object combination
				cmb2 = tricomb[i*60 + 1 * *nTri + j]; // user object combination
				cmb3 = tricomb[i*60 + 2 * *nTri + j]; // user object combination							
				d1 = x[cmb1-1];
				d2 = x[cmb2-1];
				d3 = x[cmb3-1];																																								
				// get evolutionary model stationary points
				funsimX( &x12eq, &x13eq, &d1, &d2, &d3, &Bs );
				// trained model for inherent similarity
				funsearchTriArray( combin, &cmb1, &cmb2, &cmb3, &f1r, &f1c, &f2r, &f2c, &f3r, &f3c, nTri );
				// total prob of pair 1,2
				pr12 = upar * x12eq + ( 1 - upar ) * inhersim[f3r + f3c*20]; 
				// prob of pair 1,3
				pr13 = upar * x13eq + ( 1 - upar ) * inhersim[f2r + f2c*20];
				// prob of pair 2,3
				pr23 = upar * (1 - x12eq - x13eq ) + ( 1 - upar ) * inhersim[f1r + f1c*20];
				// final log probabilities
				// first and second object grouped
				lpr12f = log( pr12 );			
				// first and third object grouped
				lpr13f = log( pr13 );			
				// second and third object grouped
				lpr23f = log( pr23 );
				
				if( y[i*60 + 0 * *nTri + j]==1 ){ 
					*candLL = *candLL + lpr23f;	
					}
				if( y[i*60 + 1 * *nTri + j]==1 ){	
					*candLL = *candLL + lpr13f;
					 }
				if( y[i*60 + 2 * *nTri + j]==1 ){ 
					*candLL = *candLL + lpr12f;
					 }
				}
			}
	"   
LLt <- signature( x="numeric", y="integer", alphacomb="integer", nObj="integer", combin="integer", tricomb="integer", inhersim="numeric", candLL="numeric", N="integer", nTri="integer", debug="numeric" )

					
TriArrayFn <- "
		for( int i=0; i<*nTri; i++ ){
				*f1r = -1;	
				*f1c = -1;
				*f2r = -1;	
				*f2c = -1;
				*f3r = -1;	
				*f3c = -1;			
			for( int j=0; j<3; j++ ){
				if( combin[j * *nTri + i] == *cmb1 ){
					*f1r = i;
					*f1c = j;						
					}
				if( combin[j * *nTri + i]	== *cmb2 ){
					*f2r = i;
					*f2c = j;						
					}
				if( combin[j * *nTri + i] == *cmb3 ){
					*f3r = i;
					*f3c = j;						
					}					
				}
			if( *f1r>-1 && *f2r>-1 && *f3r>-1 ){					
				break;
				}
			}
"
calArraySearch <- signature( combin="integer", cmb1 ="integer", cmb2 ="integer", cmb3="integer",f1r="integer",f1c="integer", f2r="integer", f2c="integer", f3r="integer",f3c="integer",nTri="integer")

# test likelihood function

fitfns <- cfunction( list( funsimX = calevolSimL, funsearchTriArray=calArraySearch, funLL = LLt ), list( evolSimL, TriArrayFn, LLtriad ), language="C", convention=".C" )

funLL <- fitfns[["funLL"]]
