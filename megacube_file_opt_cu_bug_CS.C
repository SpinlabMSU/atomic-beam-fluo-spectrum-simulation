// requires ROOT: http://root.cern.ch/
//
//    usage:            compile using ROOT by: ".L .C+" (note the plus sign)
//                                             ".x .C+"
//      By:             Jaideep Singh
//
//      Date:           2023-06-24 - quadruple checked normalization for absorption cross section and oven angular distribution, 
//      			   - added unpolarized power broadening from Steck
//      			   - added transit time broadening from Shimoda
//
//      descp:      	simulate Doppler Broadened lineshape with the possibility of power broadening
// 
//      reference:      todo: 	input file - 2022-02-22
//      			output file - 2022-02-21
//      			full output info file - 2022-02-23
//      			oven model (big ugly equation) - 2023-06-03
//      			TODO: detector integration
//      			laser model - 2022-02-27
//      			general set of atomic lines - 2022-02-22
//      			optimize loop - 2202-02-22, for baseline parameters: 6797 seconds to 685 seconds! (just by moving repetitive stuff out of the inner loop)
//				converted all while loops into for loops and shifted values to bin centers - 2023-05-30
//				optimized integration limits for laser frequency profile and Maxwell-Boltzmann profile - 2023-05-31
//				print information about laser beam size - 2023-05-31
//				added projections for completion time - 2023-05-31
//				choose locations within the frequency and velocity bins to minimize integration inaccuracies by using average value of integrand within a bin
//					fMB - 2023-05-31
//					fG - 2023-05-31
//				frequencies chosen to be finely spaced around peaks and then coarsely spaced between peaks - output needs to be sorted - 2023-05-31
//				oven nozzle exit point randomization - 2023-06-04
//					needs the nozzle_radius in meters as input
//					implemented as a Monte Carlo: each microcube has a different oven nozzle exit point that is chosen randomly uniformly with a circle
//					need to increase the integration limits and bins of the megacube to account for the spatial extent of the exit nozzle
//						the x direction is the laser direction, this is fine
//						the y direction is the detector direction, this needs to be adjusted
//				normalized the oven angular distrbution to Knorm = int(j(theta)*sin(theta)*2.0,theta = 0 to 90 deg) - 2023-06-07
//				checked angular distribution normalization - added a missing factor of pi: j(theta)/Knorm -> j(theta)/pi/Knorm - 2023-06-23
//				changed normalization for the integrated cross section using Axener 2004 notation - 2023-06-19
//					reference:   2023-06-19: calculated angular factor (Aqp = $A'_q$) using Axener 2004 notation: Spectrochimica Acta Part B 59 (2004) 1–39
//      				https://doi.org/10.1016/j.sab.2003.10.002
//      				scalexs correctly has factor of 3.0 according to Corney and Steck
//      				must use check_Aq_pol_gen.C to get FFpsum now!
//				power broadening option using Steck formulation - 2023-06-19
//					does not account for splitting of the line (Mollow triplets)
//				transit time broadening using Shimoda - 2023-06-24

// head files needed to compile macro
#include "TDatime.h"
#include "TComplex.h"
#include "TMath.h"
#include "TString.h"
#include "TRegexp.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TF1.h"
#include "TGaxis.h"
#include "TH2D.h"
#include "TLatex.h"
#include "TLine.h"
#include "TStyle.h"
#include "Riostream.h"
#include "TStopwatch.h"
#include "TVector3.h"
#include "TROOT.h"
#include "TRandom3.h"
#include "TGraphErrors.h"
#include "Math/SpecFuncMathMore.h"


// BUG START: 
// translated from Big Ugly Equation (BUG) Python - Scoles Eqns. 4.29 to 4.34	 
// gamma = nozzle diamter to nozzle length ratio = Scoles beta = aspect ratio
//
Double_t fzeta0(Double_t gamma) {
	Double_t gamma2 = gamma*gamma;
	Double_t gamma3 = gamma2*gamma;
	Double_t s1g2 = TMath::Sqrt(1.0 + gamma2);
	Double_t ash = TMath::ASinH(1.0/gamma);
	return 0.5 - (1.0/3.0/gamma2)*(1.0-2.0*gamma3+(2.0*gamma2-1.0)*s1g2)/(s1g2-gamma2*ash);
}

Double_t integrandS(Double_t z, Double_t p, Double_t deltaprime, Double_t xi1, Double_t xi0) {
    Double_t xx = deltaprime*(1.0+z/p*(xi1/xi0-1.0));
    Double_t erfxx = TMath::Erf(xx);
    Double_t erfdp = TMath::Erf(deltaprime);
    return TMath::Sqrt(1.0-z*z)*(erfxx - erfdp);
}

// Simpson's Rule Numerical Integration
// Npts = number of points, Npts - 1 = must be even and equal to twice the number of parabolic segments
Double_t fSimpson(Int_t Npts, Double_t h, Double_t *y) {
        Double_t sum = 0.0;
        for (Int_t k = 0 ; k < Npts ; k++) {
                if (k == 0 or k == Npts-1) {
                        sum = sum + y[k];
                } else if (k % 2 == 1) {
                        sum = sum + 4.0*y[k];
                } else  {
                        sum = sum + 2.0*y[k];
                }
        }
	if (TMath::IsNaN(sum*h) == 1) {
		cout << "sum = " << sum << endl;
		cout << "h = " << h << endl;
	}
        return sum*h/3.0;
}

Double_t Sintegral(Double_t zmax, Double_t p, Double_t deltaprime, Double_t xi1, Double_t xi0) {
	Int_t Npts = 10001;
	Double_t zmin = 0.0;
	Double_t dz = (zmax-zmin)/(Npts-1.0);
	Double_t y[Npts];
	for (Int_t k = 0 ; k < Npts ; k++) {
		Double_t z = zmin + dz*k;
		y[k] = integrandS(z, p, deltaprime, xi1, xi0);
		if (TMath::IsNaN(y[k]) == 1) {
			cout << z << "\t" << p << "\t" << deltaprime << "\t" << xi1 << "\t" << xi0 << endl;
		}
	}
    	return fSimpson(Npts,dz,y);
}

// Big Ugly Equation, Scoles Eqns. 4.29 to 4.34
// angular distribution function 
// theta = angle in radians
// Knl = Knudsen number, infinity for molecular flow
// gamma = nozzle diamter to nozzle length ratio = Scoles beta = aspect ratio
//
Double_t j(Double_t theta, Double_t Knl, Double_t gamma) {
	Double_t zeta0 = fzeta0(gamma);
    	Double_t zeta1 = 1.0-zeta0;
//    
    	Double_t xi0 = zeta0;
    	Double_t xi1 = zeta1;
//    
    	Double_t delta = xi0/TMath::Sqrt(2.0*Knl*(xi1-xi0));
    	Double_t deltaprime = delta/TMath::Sqrt(TMath::Cos(theta));
//    
    	Double_t xi = (xi1/xi0)*delta;
    	Double_t xiprime = (xi1/xi0)*deltaprime;
//    
    	Double_t p = TMath::Tan(theta)/gamma;
//    
	Double_t output = 0.0;
	Double_t pi = TMath::Pi();
	Double_t delta2 = delta*delta;
    	Double_t erfxi = TMath::Erf(xi);
    	Double_t erfd = TMath::Erf(delta);
	Double_t xi2 = xi*xi;
	Double_t xiprime2 = xiprime*xiprime;
	Double_t deltaprime2 = deltaprime*deltaprime;
    	Double_t erfxip = TMath::Erf(xiprime);
    	Double_t erfdp = TMath::Erf(deltaprime);

    	if (p == 0.0) {
        	output = (TMath::Sqrt(pi)/2.0)*(TMath::Exp(delta2)/delta)*(erfxi-erfd);
        	output = xi0*(1.0+output);
        	output = output + (1.0-xi1)*TMath::Exp(-xi2+delta2);
		if (TMath::IsNaN(output) == 1) {
			cout << "p = 0: " << theta << "\t" << p << "\t output = " << output << endl;
		}
	} else if (p <= 1.0) {
        	Double_t F = (2.0/TMath::Sqrt(pi))*xiprime*(1.0/xi1-1.0)*TMath::Exp(-xiprime2); ///HERE
        	Double_t Rp = TMath::ACos(p) - p*TMath::Sqrt(1.0-p*p);
        	Double_t Spp = Sintegral(p,p,deltaprime,xi1,xi0);
        	output = (2.0/TMath::Sqrt(pi))*(TMath::Exp(deltaprime2)/deltaprime)*(Rp/2.0*(erfxip-erfdp+F) + Spp);
        	output = xi0*TMath::Cos(theta)*(1.0 + output);
		if (TMath::IsNaN(output) == 1) {
			cout << "p <= 1: " << theta << "\t" << p << "\t" << F << "\t" << Rp << "\t" << Spp << endl;
		}
	} else if (p > 1.0) { 
        	Double_t S1p = Sintegral(1.0,p,deltaprime,xi1,xi0);
        	output = xi0*TMath::Cos(theta)*(1.0 + (2.0/TMath::Sqrt(pi))*(TMath::Exp(deltaprime2)/deltaprime)*S1p);
		if (TMath::IsNaN(output) == 1) {
			cout << "p > 1: " << theta << "\t" << p << "\t" << S1p << endl;
		}
	} else {
        	output = 1.0E20;
	}
    	return output;
}

// theta = angle in radians
// gamma = nozzle diamter to nozzle length ratio = Scoles beta = aspect ratio
// Knl = Knudsen number, infinity for molecular flow
Double_t integrandK(Double_t theta, Double_t Knl, Double_t gamma) {
   	return 2.0*TMath::Sin(theta)*j(theta,Knl,gamma);
}

// transmission function K = olander/Kruger, Scoles W
// Knl = Knudsen number, infinity for molecular flow
// gamma = nozzle diamter to nozzle length ratio = Scoles beta = aspect ratio
Double_t K(Double_t Knl, Double_t gamma) {
	Double_t pi = TMath::Pi();
	Double_t thetamin = 0.0;
    	Double_t thetamax = 89.9/180.0*pi;
	Int_t Npts = 10001;
	Double_t dtheta = (thetamax-thetamin)/(Npts-1.0);
	Double_t y[Npts];
	for (Int_t k = 0 ; k < Npts ; k++) {
		Double_t theta = thetamin + dtheta*k;
		y[k] = integrandK(theta, Knl, gamma);
	}
    	return fSimpson(Npts,dtheta,y);
}

// Scoles' 4.16, beta = gamma = nozzle aspect ratio
Double_t falpha(Double_t beta) {
	Double_t beta2 = beta*beta;
	Double_t beta3 = beta*beta2;
	Double_t s1b2 = TMath::Sqrt(1.0 + beta2);
	Double_t ash = TMath::ASinH(1.0/beta);
    	Double_t output = 1.0 - 2.0*(beta3) + (2.0*beta2-1.0)*s1b2;
    	output = output/(s1b2 - (beta2)*ash);
    	output = 0.5 - output/3.0/beta2;
    	return output;
}

// transmission function K = olander/Kruger, Scoles W
// gamma = nozzle diamter to nozzle length ratio = Scoles beta = aspect ratio
// Knl = Knudsen number, infinity for molecular flow
Double_t fW(Double_t beta) {
	Double_t alpha = falpha(beta);
	Double_t beta2 = beta*beta;
	Double_t s1b2 = TMath::Sqrt(1.0 + beta2);
	Double_t output = 1.0 + 2.0/3.0*(1.0-2.0*alpha)*(beta-s1b2)+2.0/3.0*(1.0+alpha)/(beta2)*(1.0-s1b2);
	return output;
}
// Python
// BUG END

// angular distribution out of the oven - unitless
//Double_t j( Double_t theta) {
//	return TMath::Cos(theta);
//}
Double_t fj( Double_t theta, Double_t *par) {
	Double_t Knl   = par[0]; // Knudsen number, >10 is molecular flow regime, <0.01 is continuum flow regime
	Double_t gamma = par[1]; // diameter to length ratio of nozzle
	return j(theta,Knl,gamma);
}

// Voigt profile - ROOT definition is fast and normalized
// sigma = Gaussian width (standard deviation definition)
// FWHM = Lorentzian FWHM
// nu0 = center frequency
Double_t fV(Double_t nu , Double_t *par) {
	Double_t sigma = par[0];
	Double_t FWHM = par[1];
	Double_t nu0 = par[2];
	return TMath::Voigt( (nu-nu0) , sigma , FWHM , 4 );
}

// Lorentzian - integral normalized to 1 - units of 1/FWHM
// FWHM = full width half maximum
// nu0 = atomic resonance frequency
// nu = laser frequency
// pbopt = power broadening option, 0 = off, 1 = on
// Ir = local laser intensity
// Isat = unpolarized saturation intensity - Steck
Double_t fL(Double_t nu, Double_t *par) {
	Double_t FWHM = par[0];
	Double_t nu0 = par[1];
	Double_t pbopt = par[2];
	Double_t Ir = par[3];
	Double_t Isat = par[4];
	Double_t pi = TMath::Pi();
	return FWHM/2.0/pi/(FWHM*FWHM/4.0*(1.0+pbopt*Ir/Isat) + (nu-nu0)*(nu-nu0));
}

// Mawell-Boltzmann - integral normalized to 1 - units of s/m
// MkT = s^2/m^2 in SI
// v = m/s
Double_t fMB(Double_t v, Double_t MkT) {
	Double_t output = TMath::Sqrt(2.0/TMath::Pi());
	output = output*TMath::Power(MkT,3.0/2.0);
	output = output*v*v*TMath::Exp(-v*v*MkT/2.0);
	return output;
}

Double_t fMBf(Double_t *x, Double_t *par) {
	Double_t v = x[0];
	Double_t MkT = par[0];
	return fMB(v,MkT);
}

// Gaussian - integral normalized to 1 - units of 1/FWHM
// FWHM = full width half maximum
// nu0 = resonance frequency
// nu = laser frequency
Double_t fG(Double_t nu, Double_t *par) {
	Double_t FWHM = par[0];
	Double_t nu0 = par[1];
	Double_t pi = TMath::Pi();
	Double_t log2 = TMath::Log(2.0);
	Double_t output = 2.0*TMath::Sqrt(log2/pi)/FWHM;
	output = output*TMath::Exp(-4.0*log2*(nu-nu0)*(nu-nu0)/FWHM/FWHM);
	return output;
}

Double_t fGf(Double_t *x, Double_t *par) {
	Double_t nu = x[0];
	return fG(nu,par);
}

// getting velocity bin locations to improve accuracy of Doppler profile integration
void get_varrayopt(Int_t Nvel , Double_t dv , Double_t MkT , Double_t *varray , Double_t *varrayopt) {
	Double_t vmin = varray[0];
	Double_t vmax = varray[Nvel-1];
	TF1 *func = new TF1("func",fMBf,vmin,vmax,1);
	Double_t par[1];
	par[0] = MkT;
	func->SetParameters(par);
	for (Int_t n = 0 ; n < Nvel ; n++) {
		vmin = varray[n] - dv/2.0;
		vmax = varray[n] + dv/2.0;
		Double_t meanfMB = func->Integral(vmin,vmax)/(vmax-vmin);
		varrayopt[n] = func->GetX(meanfMB,vmin,vmax);
	}
}

// getting nugamma bin centers based on locations of transitions
// fine spacing within +/- FWHMDB/2.0 of peak center
// coarse spacing elsewhere = (coarse)*(fine spacing)
// nu_i = array of tranistion locations, must be in numerical order for this to work
// nugarray = equally spaced frequency bins
// output array is nugtransarray = is not sorted, needs to be sorted
void get_nugtransarray(Int_t coarse , Int_t numtrans , Double_t *nu_i , Double_t FWHMDB , Int_t Ngam , Double_t *nugarray, Double_t *Ngamtrans , Double_t *nugtransarray) {
	Double_t dnu = nugarray[1] - nugarray[0];
	Double_t numin = nugarray[0];
	Double_t numax = nugarray[Ngam-1];
	Int_t k = 0;
	for(Int_t n = 0 ; n < numtrans ; n++) {
		Int_t l = 0;
		Double_t thisnu = numin;
		if (n == 0) {
			Double_t firstnumin = nu_i[n] - FWHMDB/2.0;
			while (thisnu < firstnumin) {
				nugtransarray[k] = thisnu;
				k++;
				l++;
				thisnu = numin + l*dnu*coarse;
			}
		}
		numin = nu_i[n] - FWHMDB/2.0;
		numax = nu_i[n] + FWHMDB/2.0;
		l = 0;
		thisnu = numin + l*dnu;
		while (thisnu < numax) {
			nugtransarray[k] = thisnu;
			k++;
			l++;
			thisnu = numin + l*dnu;
		}
		if (n < numtrans-1) {
			Double_t nextnumin = nu_i[n+1] - FWHMDB/2.0;
			if (nextnumin > numax) {
				l = 0;
				thisnu = numax + dnu*coarse*l;
				while(thisnu < nextnumin) {
					nugtransarray[k] = thisnu;
					k++;
					l++;
					thisnu = numax + dnu*coarse*l;
				}	
			}
		}
		if (n == numtrans - 1) {
			Double_t lastnumax = nugarray[Ngam-1];
			l = 0;
			thisnu = numax + dnu*coarse*l;
			while(thisnu < lastnumax) {
				nugtransarray[k] = thisnu;
				k++;
				l++;
				thisnu = numax + dnu*coarse*l;
			}
		}
	}
	Ngamtrans[0] = k;
}

// getting relative frequency bin locations to improve accuracy of laser profile integration
// relnumin = relative min frequency
// relnumax = relative max frequency
// gpar0 = FWHM of Gaussian
// Nnu = number of nu bins
// relnuarray = relative nu bin locations - this will be used in the freq loop after providing relative shifts to nugamma
// fGarray = value of Gaussian at relative bin locations - this will be used in the freq loop
void get_fGarrays(Double_t relnumin , Double_t relnumax , Double_t gpar0 , Int_t Nnu , Double_t *relnuarray, Double_t *fGarray) {
	Double_t dnu = (relnumax - relnumin)/Nnu;
	TF1 *func = new TF1("func",fGf,relnumin,relnumax,2);
	Double_t par[2];
	par[0] = gpar0;
	par[1] = 0.0;
	func->SetParameters(par);
	for (Int_t n = 0 ; n < Nnu ; n++) {
		Double_t numin = relnumin + n*dnu;
		Double_t numax = numin + dnu;
		Double_t meanfG = func->Integral(numin,numax)/(numax-numin);
		relnuarray[n] = func->GetX(meanfG,numin,numax);
		fGarray[n] = func->Eval(relnuarray[n]);
	}
}


// main program
// loops over laser frequency center
// 	laser spectral profile
// 		maxwell-boltzmann velocity profile
// 			x-coordinate = laser direction
// 				y-coordinate = detector direction
// 					z-coordinate = oven direction
// all units SI except for laser frequency which is MHz
// Nxyz = number of microcubes along a single directions, Nxyz***3 = number of microcubes in MEGACUBE
// bquick = true if only 11 laser frequencies are quickly scanned through, otherwise will do the whole spectrum
//
void megacube_file_opt_cu_bug_CS(TString inputfilebase = "RbD1") {
	gROOT->Reset();
	gROOT->DeleteAll(); 
// starting Timer
	TStopwatch timer;
	TDatime timestamp;
	cout << "START" << endl;
	timer.Start();
	timestamp.Set();
	timestamp.Print();
// constants
	Double_t eps = 1.0E-10; 	// small number 
	Double_t c = 299792458.0;       // m/s, speed of light
	Double_t NA = 6.02214076E23;    // particles per mol, Avogadro
	Double_t kB = 1.380649E-23;     // J/K, Boltzmann
	Double_t h = 6.62607004E-34;    // J/Hz, Planck
	Double_t re = 2.8179403262E-15; // m, classical electron radius
	Double_t pi = TMath::Pi();      // unitless, just pi
// parameters from input file - default values for Rb
	Int_t Nnu = 100; 		// integer number of laser profile integral bins, 
	Int_t Nvel = 300;		// integer number of Maxwell-Boltzmann profile integral bins
	Double_t xmin = -0.5E-2;	// megacube minimum x position in meters, 
	Double_t xmax = +0.5E-2;	// megacube maximum x position in meters, 
	Double_t ymin = -0.5E-2;	// megacube minimum y position in meters, 
	Double_t ymax = +0.5E-2;	// megacube maximum y position in meters, 
	Double_t zmin = -0.5E-2;	// megacube minimum z position in meters, 
	Double_t zmax = +0.5E-2;	// megacube maximum z position in meters, 
	Int_t Nxtot = 6;		// integer number of spatial intervals along the x direction = laser direction, 
	Int_t Nytot = 6;		// integer number of spatial intervals along the y direction = detector direction, 
	Int_t Nztot = 6;		// integer number of spatial intervals along the z direction = oven direction, 
	Bool_t bovenrandom = kTRUE;	// boolean that determines whether the Doppler Broadening is completely Random as in a vapor cell (kTRUE) or based on geometry (kFALSE), 
	Int_t ranseed = 23;		// integer seed for random variable, used if kTRUE above, 
	Bool_t bquick = kTRUE;		// boolean which indicates whether the signal is calculated for just a few laser frequencies (kTRUE) or for the full set of frequencies (kFALSE), 
	Double_t fa = 0.342;		// oscillator strength (unitless), 
	Double_t Anat = 3.61E7;		// transition natural linewidth in rad*Hz, Einstein A, 
					//
					//
					// this is where the resonant frequency array stuff is entered
					//
					//
	Double_t LP = 1.0E-4;		    // laser power in Watts, 
	Double_t waist = 1.0E-2;	    // laser beam waist in meters, 
	Double_t waist_locationx = 0.0;	    // laser beam waist location along x direction in meters, 
	Double_t laser_linewidth = 5.0;	    //  laser linewidth FWHM in MHz, 
	Double_t laser_step = 4.0/5.0;	    //  laser step size in units of laser linewidth, unitless 
	TVector3 detloc(+0.0,+0.0774,+0.0); //  detector location in meters, tab-separated vector (x,y,z), 
	Double_t detx = 1.0E-3; 	    // detector dimensions (assuming rectangle) in meters, tab-separated vector (x,z), 
	Double_t detz = 1.0E-3;		    // detector dimensions (assuming rectangle) in meters, tab-separated vector (x,z), 
	TVector3 ovnloc(+0.0,+0.0,-0.132);  //  oven location in meters, tab-separated vector (x,y,z), 
	Double_t T = 100.0+273.15;          // Kelvin, oven temperature
	Double_t Kn = 1.0E10;		    //  oven Knudsen number, unitless, 
	Double_t gamma_oven = 1.0E10;       //  oven nozzle diameter to to length ratio, unitless, 
	Double_t nozzle_radius = 0.0016;    // nozzle exit radius, meters, default = 0.125 inch diameter (short Yb nozzle)
	Double_t pbopt = 0.0;		    // power broadening option, 0 = off, 1 = on
	Double_t ttopt = 0.0;		    // transit time broadening option, 0 = off, 1 = on
// reading in parameters from input file
	TString inputfilename = inputfilebase + ".txt";
	ifstream inputfile;
	inputfile.open(inputfilename);
	inputfile >> Nnu;
	cout << "Nnu = " << Nnu << endl;
	inputfile >> Nvel;
	cout << "Nvel = " << Nvel << endl;
	inputfile >> xmin;
	cout << "xmin = " << xmin << endl;
	inputfile >> xmax;
	cout << "xmax = " << xmax << endl;
	inputfile >> ymin;
	cout << "ymin = " << ymin << endl;
	inputfile >> ymax;
	cout << "ymax = " << ymax << endl;
	inputfile >> zmin;
	cout << "zmin = " << zmin << endl;
	inputfile >> zmax;
	cout << "zmax = " << xmax << endl;
	inputfile >> Nxtot;
	cout << "Nxtot = " << Nxtot << endl;
	inputfile >> Nytot;
	cout << "Nytot = " << Nytot << endl;
	inputfile >> Nztot;
	cout << "Nztot = " << Nztot << endl; 
	inputfile >> bovenrandom;
	cout << "bovenrandom = " << bovenrandom << endl;
	inputfile >> ranseed;
	cout << "ranseed = " << ranseed << endl;
	inputfile >> bquick;
	cout << "bquick = " << bquick << endl;
	inputfile >> fa;
	cout << "fa = " << fa << endl;
	inputfile >> Anat;
	cout << "Anat = " << Anat << endl;
	//
	//
	// this is where the resonant frequency array stuff is entered
	//
	//
	Int_t in_numtrans = 8;
	inputfile >> in_numtrans;
	cout << "num freq = " << in_numtrans << endl;
	const Int_t numtrans = in_numtrans; 	// integer number of resonant frequecies rows (see below), 
	cout << "num freq check = " << numtrans << endl;
	Double_t M_i[numtrans];			// array of masses (g/mol),
	Double_t nu_i[numtrans];		// array of absolute frequencies (MHz),
	Double_t w_i[numtrans]; 		// array of weight factors (unitless)
// weight factor = (eta_i)*(sum_{mF} rho_{F,mF})*(p_q) 
// eta_i = fraction of isotope i, sum_i eta_i = 1, natural abundance is assumed
// rho_{F,mF} = relative population of ground state, sum_F sum_{mF} rho_{F,mF} = 1, equal ground state population is assumed =  1/(2J+1)/(2I_i+1) where J is ground total electronic and I_i is the nuclear spin of isotope i
// 	in this case sum_{mF} rho_{F,mF} = (2F+1)/(2J+1)/(2I_i+1) 
// p_q = magnitude of light polarization component q, sum_q p_q = 1, p_q = 1/3 for unpolarized light which is assumed
// weights are generated by multipying FFPsum (from check_Aq_pol_gen.C) by eta_i 
	for (Int_t i = 0 ; i < numtrans ; i++) {
		inputfile >> M_i[i] >> nu_i[i] >> w_i[i];
	}
	for (Int_t i = 0 ; i < numtrans ; i++) {
		cout << "frequency #" << i+1 << "\t" <<  M_i[i] << "\t" << nu_i[i] << "\t" <<  w_i[i] << endl;
	}
	//
	//
	// this is where the resonant frequency array stuff is entered
	//
	//
	inputfile >> LP;
	cout << "LP = " << LP << endl;
	inputfile >> waist;
	cout << "waist = " << waist << endl;
	inputfile >> waist_locationx;
	cout << "waist_locationx = " << waist_locationx << endl;
	inputfile >> laser_linewidth;
	cout << "laser_linewidth = " << laser_linewidth << endl;
	inputfile >> laser_step;
	cout << "laser_step " << laser_step <<  endl;
	Double_t x,y,z;						// generic position variables
	inputfile >> x >> y >> z;
	detloc.SetXYZ(x,y,z);
	detloc.Print();
	inputfile >> detx >> detz;
	cout << "detx = " << detx << endl;
	cout << "detz = " << detz << endl;
	inputfile >> x >> y >> z;
	ovnloc.SetXYZ(x,y,z);
	ovnloc.Print();
	inputfile >> T;
	cout << "T = " << T << endl;
	inputfile >> Kn;
	cout << "Kn = " << Kn << endl;
	inputfile >> gamma_oven;
	cout << "gamma_oven = " << gamma_oven << endl;
	inputfile >> nozzle_radius;
	cout << "nozzle radius = " << nozzle_radius << endl;
	inputfile >> pbopt;
	cout << "power broadening = " << pbopt << endl;
	inputfile >> ttopt;
	cout << "transit time broadening = " << ttopt << endl;
	inputfile.close();
// velocity factors for Maxwell-Boltzmann distribution
	Double_t MkT_i[numtrans];  	// array of velocity factors for each transition
	Double_t MkT = 1.0;             // geometric mean of velocity factors
	for (Int_t i = 0 ; i < numtrans ; i++) {
		MkT_i[i] = M_i[i]/1000.0/kB/T/NA;
		MkT = MkT*MkT_i[i];
	}
	MkT = TMath::Power(MkT,1.0/numtrans);                   // (m/s)^2
// setting the integral limits for the MB distribution integral 
	Double_t FWHMDB = TMath::Sqrt(8.0/MkT*TMath::Log(2.0)); // m/s
	Double_t vmax = TMath::Sqrt(2.0/MkT) + 1.5*FWHMDB;      // m/s
	Double_t vmin = vmax/1.0E4; 				// m/s
	Double_t vmean = TMath::Sqrt(8.0/MkT/pi);		// m/s
	Double_t dv = (vmax-vmin)/Nvel;				// m/s
	Double_t varray[Nvel],varrayopt[Nvel];
	for (Int_t l = 0 ; l < Nvel ; l++) {
		varray[l] = vmin + dv/2.0 + l*dv;		// m/s
	}
        get_varrayopt(Nvel , dv , MkT , varray , varrayopt);
// checking weight sum
	Double_t wsum = 0.0;   // weight check sum
	for (Int_t i = 0 ; i < numtrans ; i++) {
		wsum += w_i[i];					// unitless
	}
// laser parameters - frequencies in MHz
// setting the laser scan range based on the maximum Doppler broadening
	Double_t gpar[2];
	gpar[0] = laser_linewidth;      // laser linewidth MHz
	Double_t nugamma0 = 0.0;	// average transition frequency, MHz
	for (Int_t i = 0 ; i < numtrans ; i++) {
		nugamma0 += w_i[i]*nu_i[i]/wsum;
	}
	FWHMDB = FWHMDB*nugamma0/c; // MHz
	Double_t nugammamin = TMath::MinElement(numtrans,nu_i)-1.5*FWHMDB; // MHz
	Double_t nugammamax = TMath::MaxElement(numtrans,nu_i)+1.5*FWHMDB; // MHz
	Double_t dnugamma = gpar[0]*laser_step;  // MHz
	Int_t coarse = 10; // ratio of coarse spacing to fine spacing
// setting up the laser frequency scan array
	Int_t knug = 0;
	const Int_t Ngam = TMath::Nint((nugammamax-nugammamin)/dnugamma) + 1;
	Double_t nugarray[Ngam];
	for (knug = 0 ; knug < Ngam; knug++) {
		nugarray[knug] = nugammamin + knug*dnugamma;
	}
	Int_t Nmax = Ngam;
	Double_t Ngamtrans[1];
	Double_t nugtransarray[Ngam*2];
// quick scan are the peaks and the max, mid, min frequecies
	if (bquick == kTRUE) {
		nugarray[0] = nugammamin;
		nugarray[1] = nugammamax;
		nugarray[2] = nugamma0;
		for (knug = 3 ; knug < 3 + numtrans ; knug++) {
			nugarray[knug] = nu_i[knug-3];
		}
		Nmax = 3 + numtrans;
		for (Int_t k = 0 ; k < Nmax ; k++) {
			nugtransarray[k] = nugarray[k];
		}
	} else {
		get_nugtransarray(coarse , numtrans , nu_i , FWHMDB , Ngam , nugarray, Ngamtrans , nugtransarray);
		Nmax = TMath::Nint(Ngamtrans[0]);
	}
// atomic parameters
	Double_t Isat = Anat*Anat*h*(nugamma0*1.0E6)/8.0/pi/re/c/fa;  // Hz*Hz*J*s*MHz*Hz/MHz/m/(m/s), W/m^2, unpolarized saturation intensity 
	cout << "Isat = " << Isat << endl;
	Double_t lpar[5];
	lpar[0] = Anat/2.0/pi/1.0E6;  // Einstein A for Rb, MHz
	lpar[2] = pbopt; // power broadening option, 0 = off, 1 = on
	lpar[4] = Isat; // unpolarized saturation intensity W/m^2
	Double_t vpar[3];
// megacube parameters
	Double_t dx = (xmax - xmin)/Nxtot; // m
	Double_t dy = (ymax - ymin)/Nytot; // m
	Double_t dz = (zmax - zmin)/Nztot; // m
// location parameters
	TVector3 xlaser(1.0,0.0,0.0); // laser axis, unit vector, unitless
	TVector3 ydet(0.0,1.0,0.0); // detector axis, unit vector, unitless
	TVector3 zoven(0.0,0.0,1.0); // oven axis, unit vector, unitless
// bovenrandom selects the angle of the atom velocity vector randomly - True or False = wrt to oven to microcube displacement vector 	
	TRandom3 rnum(ranseed);
// Setting up output files
	TString firstname = "DB";
	if (pbopt > 0.0) {
		firstname = firstname + "PB";
	}
	if (ttopt > 0.0) {
		firstname = firstname + "TT";
	}
	firstname = firstname + "-";
	TString basename = TString::Format("-%03dMC-%03d-%03d-%03d-%03d-%03d-%03d",TMath::Nint((xmax-xmin)*1000),TMath::Nint(dnugamma),Nnu,Nvel,Nxtot,Nytot,Nztot);
	basename = firstname + inputfilebase + basename;
	if (bovenrandom == kTRUE) {
		basename += "-Random";
	} else {
		basename += "-BUG";
	}
	if (bquick == kTRUE) {
		basename += "-quick";
	}
	TString outputfilename = basename + ".txt";
	TString outinfofilename = basename + "-info.txt";
	ofstream outputfile;
	ofstream outinfo;
	outputfile.open(outputfilename);
	outinfo.open(outinfofilename);
	cout << scientific;
	cout << setprecision(11);
	outputfile << scientific;
	outputfile << setprecision(11);
	outinfo << scientific;
	outinfo << setprecision(11);
	cout << outputfilename << endl;
	cout << outinfofilename << endl;
	outinfo << outputfilename << endl;
	outinfo << outinfofilename << endl;
	outinfo << "Nnu = " << Nnu << endl;
	outinfo << "Nvel = " << Nvel << endl;
	outinfo << "xmin = " << xmin << endl;
	outinfo << "xmax = " << xmax << endl;
	outinfo << "ymin = " << ymin << endl;
	outinfo << "ymax = " << ymax << endl;
	outinfo << "zmin = " << zmin << endl;
	outinfo << "zmax = " << xmax << endl;
	outinfo << "Nxtot = " << Nxtot << endl;
	outinfo << "Nytot = " << Nytot << endl;
	outinfo << "Nztot = " << Nztot << endl;
	outinfo << "bovenrandom = " << bovenrandom << endl;
	outinfo << "ranseed = " << ranseed << endl;
	outinfo << "bquick = " << bquick << endl;
	outinfo << "fa = " << fa << endl;
	outinfo << "Anat = " << Anat << endl;
	for (Int_t i = 0 ; i < numtrans ; i++) {
		cout    << "vmp (" << i << ") = " << TMath::Sqrt(2.0/MkT_i[i]) << " m/s" << endl;
		outinfo << "vmp (" << i << ") = " << TMath::Sqrt(2.0/MkT_i[i]) << " m/s \t frequency #" << i+1 << "\t" <<  M_i[i] << "\t" << nu_i[i] << "\t" <<  w_i[i] << endl;
	}
	cout << "v range = " << vmin << "\t" << vmean << "\t" << vmax << " m/s" << endl;
	cout << "wsum = " << wsum << " (should be 1/3 for unpolarized light)" << endl;
	cout << "lambda0 = " << c/(nugamma0*1E6)/1E-9 << " nm" << endl;
	cout << "DB FWHM = " << FWHMDB << " MHz" << endl;
	cout << "nu gamma min detuning = " << nugammamin - nugamma0 << " MHz" << endl;
	cout << "nu gamma max detuning = " << nugammamax - nugamma0 << " MHz" << endl;
	outinfo << "v range = " << vmin << "\t" << vmean << "\t" << vmax << " m/s" << endl;
	outinfo << "wsum = " << wsum << " (should be 1)" << endl;
	outinfo << "lambda0 = " << c/(nugamma0*1E6)/1E-9 << " nm" << endl;
	outinfo << "DB FWHM = " << FWHMDB << " MHz" << endl;
	outinfo << "nu gamma min detuning = " << nugammamin - nugamma0 << " MHz" << endl;
	outinfo << "nu gamma max detuning = " << nugammamax - nugamma0 << " MHz" << endl;
	outinfo << "LP = " << LP << endl;
	outinfo << "waist = " << waist << endl;
	outinfo << "waist_locationx = " << waist_locationx << endl;
	outinfo << "laser_linewidth = " << laser_linewidth << endl;
	outinfo << "laser_step " << laser_step <<  endl;
	cout << "number of nu gammas = " << Nmax << endl;
	outinfo << "number of nu gammas = " << Nmax << endl;
	outinfo << "detx = " << detx << endl;
	outinfo << "detz = " << detz << endl;
	outinfo << "T = " << T << endl;
	outinfo << "Kn = " << Kn << endl;
	outinfo << "gamma_oven = " << gamma_oven << endl;
	outinfo << "nozzle radius = " << nozzle_radius << endl;
	outinfo << "power broadening = " << pbopt << endl;
	outinfo << "unpolarized saturation intensity = " << Isat << " W/m^2" << endl;
	outinfo << "transit time broadening = " << ttopt << endl;
// seting up oven parameters
	Double_t ovenpar[2];
	ovenpar[0] = Kn;
	ovenpar[1] = gamma_oven;
// setting up fMB array (velocity profile)
	Double_t checknorm = 0.0;
	Double_t fMBarray[Nvel][numtrans];
	for (Int_t l = 0 ; l < Nvel ; l++) {
		outinfo << "vel = " << l << "\t" << varray[l] << "\t" << varrayopt[l] << endl;
		for (Int_t i = 0 ; i < numtrans; i++) {
			fMBarray[l][i] = fMB(varrayopt[l],MkT_i[i]);  // s/m
			checknorm = checknorm + dv*fMBarray[l][i]; // unitless
		}
	}
	checknorm = checknorm/numtrans;
	cout << "fMB normalization check = " << checknorm << " (should be 1)" << endl;
	outinfo << "fMB normalization check = " << checknorm << " (should be 1)" << endl;
//
// setting up frequency and velocity independant calculations
//
	Double_t M2 = 1.1;
	Double_t lambda0 = c/nugamma0/1.0E6; // wavelength in m
	Double_t phi0 = M2*lambda0/pi/waist;
	Double_t zR = waist/phi0;
	cout << "phi0 = " << phi0 << endl;
	outinfo << "phi0 = " << phi0 << endl;
	cout << "zR (m) = " << zR << endl;
	outinfo << "zR (m) = " << zR << endl;
	Double_t scalexs = 3.0*pi*re*c*fa/1.0E6; // cross section parameters (factor of 1/3 for unpolarized light included in input weights), m^2*Hz to m^2*MHz, moved into outerloop - correct according to Corney and Steck
//	Double_t scalexs = pi*re*c*fa/1.0E6; // cross section parameters, m^2*Hz to m^2*MHz, moved into outerloop
	Double_t Adet = detx*detz; // m^2
	Double_t jarray[Nxtot][Nytot][Nztot];
	Double_t thetagamma_array[Nxtot][Nytot][Nztot];
	Double_t Irarray[Nxtot][Nytot][Nztot];
	Double_t sigmattarray[Nxtot][Nytot][Nztot];
	Double_t scale_array[Nxtot][Nytot][Nztot];
	Double_t Volume_check = 0.0;
	Double_t Knorm = K(Kn,gamma_oven); // normalize oven angular distribution , unitless
	cout << "Knorm = int(j(theta)*sin(theta)*2.0,theta = 0 to 90 deg) = " << Knorm << endl;
	outinfo << "Knorm = int(j(theta)*sin(theta)*2.0,theta = 0 to 90 deg) = " << Knorm << endl;
	Double_t Irmaxrel = 0.0;
	Double_t sigmattmax = 0.0;
        for (Int_t Nx = 0 ; Nx < Nxtot ; Nx++) {
		x = xmin + dx/2.0 + Nx*dx;
        	for (Int_t Ny = 0 ; Ny < Nytot ; Ny++) {
			y = ymin + dy/2.0 + Ny*dy;
        		for (Int_t Nz = 0 ; Nz < Nztot ; Nz++) {
				z = zmin + dz/2.0 + Nz*dz;
				TVector3 microcube(x,y,z);
// choosing a random location within the nozzle aperature
				Double_t rnozzle = nozzle_radius*TMath::Sqrt(rnum.Uniform(1.0));
				Double_t phinozzle = 2.0*pi*rnum.Uniform(1.0);
				Double_t xnozzle = rnozzle*TMath::Cos(phinozzle);
				Double_t ynozzle = rnozzle*TMath::Sin(phinozzle);
				TVector3 thisovnloc( xnozzle + ovnloc.X() , ynozzle + ovnloc.Y() , ovnloc.Z() );
				TVector3 o2m = microcube-thisovnloc;
// calculating each laser ray from Gaussian Beam moving along x-axis
				Double_t rlaser = TMath::Sqrt(y*y+z*z);
				Double_t beam_radius = waist*TMath::Sqrt(1.0+TMath::Power(phi0*(x-waist_locationx)/waist,2.0));
				Double_t Roc = (x-waist_locationx)*(1.0 + TMath::Power(zR/(x-waist_locationx),2.0));
				sigmattarray[Nx][Ny][Nz] = (1.0/2.0/pi)*TMath::Sqrt(TMath::Power(pi*beam_radius/Roc/lambda0,2.0) + TMath::Power(beam_radius,-2.0));
				Double_t tan_ray = phi0*rlaser/beam_radius/TMath::Sqrt(1.0 + TMath::Power(waist/phi0/(x-waist_locationx+eps),2.0));
				xlaser.SetXYZ(1.0,y*tan_ray/rlaser,z*tan_ray/rlaser);
				Double_t thetagamma = o2m.Angle(xlaser);
				thetagamma_array[Nx][Ny][Nz] = thetagamma;
				Double_t thetaoven  = o2m.Angle(zoven);
				jarray[Nx][Ny][Nz] = (fj(thetaoven,ovenpar)/pi/Knorm)/o2m.Dot(o2m)*dx*dy*dz;                            // m^3/m^2 = m
				TVector3 m2d = detloc - microcube;                                                                      // m
				Double_t scale = scalexs*jarray[Nx][Ny][Nz];								// oven parameters, m^2*MHz*m = m^3*MHz
				scale = scale*Adet*ydet.Dot(m2d)/4.0/pi/TMath::Power(m2d.Mag(),3.0);					// detector parameters, m^3*MHz*m^2*m/m^3 = m^3*MHz
				Double_t Ir = LP*2.0/pi/beam_radius/beam_radius*TMath::Exp(-2.0*rlaser*rlaser/beam_radius/beam_radius);	// local laser intensity, W/m^2
				Irarray[Nx][Ny][Nz] = Ir;
				scale = scale*Ir;											// laser parameters, m^3*MHz*W/m^2 = m*MHz*W
				scale_array[Nx][Ny][Nz] = scale;                                                                        // m*MHz*W
				outinfo << x << "\t" << y << "\t" << z << "\t" << thetaoven*180.0/pi << "\t" << thetagamma*180.0/pi << "\t"; 	// x,y,z in meters ; angles in degrees
				outinfo << beam_radius << "\t" << jarray[Nx][Ny][Nz] << "\t"; 							// beam radius in meters ; j(theta)/r^2*dV in meters
			        outinfo << xnozzle << "\t" << ynozzle; 										// nozzle location in meters
			        outinfo << "\t" << Irarray[Nx][Ny][Nz]/Isat; 									// local laser intensity in units of Isat
				outinfo << "\t" << 1.0/sigmattarray[Nx][Ny][Nz] << endl;							// effective interaction length in m
				if (Ir/Isat > Irmaxrel) {
					Irmaxrel = Ir/Isat;
				}
				if ( vmean*sigmattarray[Nx][Ny][Nz] > sigmattmax ) {
					sigmattmax = vmean*sigmattarray[Nx][Ny][Nz];
				}
				Volume_check = Volume_check + dx*dy*dz;
			}
		}
	}
	cout << "max local laser intensity in units of Isat = " << Irmaxrel << endl;
	outinfo << "max local laser intensity in units of Isat = " << Irmaxrel << endl;
	cout << "max transit time linewdith using vmean = " << sigmattmax/1.0E6 << " MHz" << endl;
	outinfo << "max transit time linewdith using vmean = " << sigmattmax/1.0E6 << " MHz" << endl;
	cout << "volume = " << (xmax-xmin)*(ymax-ymin)*(zmax-zmin) << " =? " << Volume_check << " m^3" << endl; 
	outinfo << "volume = " << (xmax-xmin)*(ymax-ymin)*(zmax-zmin) << " =? " << Volume_check << " m^3" << endl; 
//
// finished setting up frequency and velocity independant calculations
//
// START laser frequency loop 	
	cout << "\n" << "########## Starting Frequency Loop\n" << endl; 
	TDatime t0;
	t0.Set();
	UInt_t t0number = t0.Convert();
	t0.Print();
	TDatime tnow;
// laser profile integral limits based on laser freq and linewidth
	Double_t relnumin = -1.5*gpar[0]; 
	Double_t relnumax = +1.5*gpar[0];
	Double_t dnu = (relnumax - relnumin)/Nnu;
// setting up fG array (laser profile)
	Double_t relnuarray[Nnu];
	Double_t fGarray[Nnu];
	get_fGarrays(relnumin , relnumax , gpar[0] , Nnu , relnuarray, fGarray);
	checknorm = 0.0;
	for (Int_t k = 0 ; k < Nnu ; k++) {
		checknorm = checknorm + dnu*fGarray[k];
		outinfo << "freq = " << k << "\t" << relnuarray[k] << "\t" << fGarray[k] << endl;
	}
	cout << "fG normalization check for nugamma = " << checknorm << " (should be 1)" << endl;
	outinfo << "fG normalization check for nugamma = " << checknorm << " (should be 1)" << endl;
	Int_t knugout = TMath::Nint(Nmax/10.0);
	for (knug = 0 ; knug < Nmax ; knug++) {
		Double_t nugamma = nugtransarray[knug];
		gpar[1] = nugamma;
		Double_t sum = 0.0;
// START frequency integral over laser spectral profile
		for (Int_t k = 0 ; k < Nnu ; k++) {
			Double_t nu = nugamma+relnuarray[k];
			Int_t l = 0;
// START velocity integral over Maxwell-Boltzmann profile
			for (Int_t l = 0 ; l < Nvel ; l++) {
				Double_t v = varrayopt[l];
// START x-integral along laser direction
        			for (Int_t Nx = 0 ; Nx < Nxtot ; Nx++) {
// START y-integral along detector direction
        				for (Int_t Ny = 0 ; Ny < Nytot ; Ny++) {
// START z-integral along oven direction
        					for (Int_t Nz = 0 ; Nz < Nztot ; Nz++) {			// converted all function calls into arrays evaluated outside loop
							Double_t thetagamma = thetagamma_array[Nx][Ny][Nz];
							Double_t Ir = Irarray[Nx][Ny][Nz];
							lpar[3] = Ir;
							vpar[0] = v*sigmattarray[Nx][Ny][Nz]/1.0E6; // transit time linewidth in MHz
							vpar[1] = Anat/2.0/pi/1.0E6*TMath::Sqrt(1.0 + pbopt*Ir/Isat); // power broadened linewidth, MHz
							Double_t nu0 = vpar[2];
							if (bovenrandom == kTRUE) {
								thetagamma = rnum.Uniform(0.0,pi);
							}
// START sum over resonant transitions			
							for (Int_t i = 0 ; i < numtrans ; i++) {
								lpar[1] = nu_i[i]*(1.0-v/c*TMath::Cos(thetagamma));
								vpar[2] = lpar[1];
								if (ttopt < 0.5) {
									sum = sum + scale_array[Nx][Ny][Nz]*w_i[i]*(dnu*fGarray[k]*fL(nu,lpar))*(fMBarray[l][i]*dv/v);  // m*MHz*W*(s/m)*(m/s)*MHz/MHz/MHz*(s/m) = m*W*(s/m) = W*s = J
								} else {
									sum = sum + scale_array[Nx][Ny][Nz]*w_i[i]*(dnu*fGarray[k]*fV(nu,vpar))*(fMBarray[l][i]*dv/v)/TMath::Sqrt(1.0+pbopt*Ir/Isat);  // m*MHz*W*(s/m)*(m/s)*MHz/MHz/MHz*(s/m) = m*W*(s/m) = W*s = J
								}
							}
						}
					}
				}
			}
		}
		outputfile << knug << "\t" << nugamma << "\t" << sum << endl; // unitless \t absolute frequency MHz \t photon energy J : note that when multiplied atoms per second, energy*(#/sec) = Watts
		if ( knug % knugout == 0 ) {
			tnow.Set();
			tnow.Print();
			UInt_t nownumber = tnow.Convert();
			TDatime tprojected;
			Double_t projectednumber = (nownumber*1.0 - t0number*1.0)/((knug+1.0)/(Nmax*1.0))+t0number*1.0;
			tprojected.Set(TMath::Nint(projectednumber));
			tprojected.Print();
			cout << TString::Format("%0.3f percent done!",TMath::Nint(100000.0*(knug+1.0)/(Nmax*1.0))/1000.0) << endl;
			outinfo << TString::Format("%0.3f percent done!",TMath::Nint(100000.0*(knug+1.0)/(Nmax*1.0))/1000.0) << endl;
		}
	}
// END laser frequency loop
// finished
	timer.Stop();
	timer.Print();
	outinfo << "Real time (sec) = " << timer.RealTime() << endl;
	outinfo << " CPU time (sec) = " << timer.CpuTime() << endl;
	outputfile.close();
	outinfo.close();
	cout << "END" << endl; 
	timestamp.Set();
	timestamp.Print();
}
