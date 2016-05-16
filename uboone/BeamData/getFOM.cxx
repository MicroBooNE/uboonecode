#include "getFOM.h"
#include <iostream>

#include <math.h>

#include "TF1.h"
#include "TRandom.h"
#include "TSpectrum.h"

using namespace gov::fnal::uboone::datatypes;

namespace bmd {
int NSTEP=1000;
int NSIGMA=3;
float TGTR=5.07;
float SIG3=0.988891003;
int NPEAKS=3;
}

float bmd::getFOM(std::string beam, const ub_BeamHeader& bh, const std::vector<ub_BeamData>& bd)
{
  float fom = -2.0;	// initialize to bogus value;
  
  if ( beam=="bnb" ) {
	bool good_tor860 = false;
	bool good_tor875 = false;
	bool good_hp875  = false;
	bool good_vp875  = false;
	bool good_hptg1  = false;
	bool good_vptg1  = false;
	bool good_mtgth  = false;
	bool good_mtgtv  = false;

	std::vector<double> tor860(1);
	std::vector<double> tor875(1);

	std::vector<double> hp875(1);
	std::vector<double> vp875(1);
	std::vector<double> hptg1(1);
	std::vector<double> vptg1(1);

	std::vector<double> mtgth(48);
	std::vector<double> mtgtv(48);


	for(auto& bdata : bd)	// get toroid, BPM, and multiwire data
	{
		if(bdata.getDeviceName().find("E:TOR860") != std::string::npos)	// get E:TOR860 reading
		{
			tor860[0] = bdata.getData()[0];
			good_tor860 = true;
			continue;
		}
		else if(bdata.getDeviceName().find("E:TOR875") != std::string::npos)	// get E:TOR875 reading
		{
			tor875[0] = bdata.getData()[0];
			good_tor875 = true;
			continue;
		}
		else if(bdata.getDeviceName().find("E:HP875") != std::string::npos)	// get E:HP875 reading
		{
			hp875[0] = bdata.getData()[0];
			good_hp875 = true;
			continue;
		}
		else if(bdata.getDeviceName().find("E:VP875") != std::string::npos)	// get E:VP875 reading
		{
			vp875[0] = bdata.getData()[0];
			good_vp875 = true;
			continue;
		}
		else if(bdata.getDeviceName().find("E:HPTG1") != std::string::npos)	// get E:HPTG1 reading
		{
			hptg1[0] = bdata.getData()[0];
			good_hptg1 = true;
			continue;
		}
		else if(bdata.getDeviceName().find("E:VPTG1") != std::string::npos)	// get E:VPTG1 reading
		{
			vptg1[0] = bdata.getData()[0];
			good_vptg1 = true;
			continue;
		}
		else if(bdata.getDeviceName().find("E:MMBTBB") != std::string::npos)	// get horizontal and vertical target multiwire reading
		{
			for( int i=0; i<48; ++i )
			{
				mtgth[i] = -bdata.getData()[i];
				mtgtv[i] = -bdata.getData()[i+48];
			}
			good_mtgth = true;
			good_mtgtv = true;
			continue;
		}
	}

	if( !(good_tor860 && good_tor875 && good_hp875 && good_vp875 && good_hptg1 && good_vptg1 && good_mtgth && good_mtgtv) )
	{
		return 2.0;	// check that all data is here; if not, return bogus value
	}

	/* add offset.  this should be a database call */
	hp875[0] -= -3.401;
	vp875[0] -=  1.475;
	hptg1[0] -=  0.457;
	vptg1[0] -=  0.389;

	/* extrapolate position to target multiwire */
	double mtgtx = hp875[0] + (hptg1[0]-hp875[0])/2.755*3.790;
	double mtgty = vp875[0] + (vptg1[0]-vp875[0])/2.348*3.586;

	/* filter multiwires */
	for( int i=0; i<48; ++i )
	{
		if( mtgth[i] < 0.0 ) mtgth[i]=0.0;
	}
	for( int i=0; i<48; ++i )
	{
		if( mtgtv[i] < 0.0 ) mtgtv[i]=0.0;
	}

	/* get profile width */
	double mtgthw = profilewidth(mtgth);
	double mtgtvw = profilewidth(mtgtv);

	/* find if edge of beam is off-target */
	int checktgt;
	double edgex;
	double edgey;
	if( mtgthw>0.0 && mtgtvw>0.0 )	// found a good width in each plane
	{
		edgex = findedge(mtgtx, mtgthw);
		edgey = findedge(mtgty, mtgtvw);

		if( edgex<TGTR && edgey<TGTR )
		{
			checktgt = 1;	// beam is within target
			fom = 1.0;
		}
		else
		{
			checktgt = 0;	// beam is not within target
		}
	}
	else
	{
		return 3.0;	// did not find a good width; return bogus value
	}

	/* this next section calculates fraction of beam off target */
	if( checktgt==0 )
	{
		int npeaks=0;
		double par[11];
		double sigmax;
		double sigmay;

		double fmt;	/* fraction missing target */

		profilesigma( mtgth, &npeaks, par );	// fit beam to gauss + linear background
		sigmax = par[4];

		profilesigma( mtgtv, &npeaks, par );
		sigmay = par[4];

		edgex = findedge(mtgtx, NSIGMA*sigmax);
		edgey = findedge(mtgty, NSIGMA*sigmay);

		if( edgex<TGTR && edgey<TGTR )	// check if beam is contained within target
		{
			fom = 1.0;
		}
		else	// if not, calculate fraction missing the target
		{
			fmt = fractionmisstarget( mtgtx, mtgty, sigmax, sigmay);
			fom = 1.0 - fmt;
		}
	}
  } else if (beam=="numi") {
    fom=100;
  }

  return fom;
}

double bmd::fractionmisstarget( double off_x, double off_y, double sig_x, double sig_y)
{
	int    i,j;	/* generic variables */

	double testR,testR2,rc2;	/* check if beam is outside target */
	double f,sum,flat;			/* summation */

	testR = TGTR;	/* radius of target */
	testR2 = testR*testR;

//	std::cout << "off_x=" << off_x << "\toff_y=" << off_y << std::endl;
//	std::cout << "sig_x=" << sig_x << "\tsig_y=" << sig_y << std::endl;

	double x,y;
#if 0
	double xmin,ymin;	/* this set is used in cartesian integral */
	double dx,dy;
	/* calculate integral using cartesian coordinates */
	xmin = bmd::NSIGMA*sig_x;
	ymin = bmd::NSIGMA*sig_y;
	dx = xmin/((double)(bmd::NSTEP));
	dy = ymin/((double)(bmd::NSTEP));
	xmin *= -1.0;
	ymin *= -1.0;

	sum=0.0;
	flat=0.0;
	x=xmin;
	for( i=0; i<2*bmd::NSTEP; ++i )
	{
		y = ymin;
		for( j=0; j<2*bmd::NSTEP; ++j )
		{
		  if( (x*x)/(bmd::NSIGMA*bmd::NSIGMA*sig_x*sig_x)+(y*y)/(bmd::NSIGMA*bmd::NSIGMA*sig_x*sig_x) < 1.0 )	/* this is three-sigma check */
			{
				rc2 = (x-off_x)*(x-off_x)+(y-off_y)*(y-off_y);
				if( rc2>testR2)
				{
					f = exp(-(((x*x)/(2.0*sig_x*sig_x))+((y*y)/(2.0*sig_y*sig_y))))/(2.0*M_PI*sig_x*sig_y);
					sum += f;
					flat += dx*dy;
				}
			}
			y += dy;
		}
		x += dx;
	}
	sum *= (dx*dy);
	flat *= (1.0-SIG3)/(M_PI*bmd::NSIGMA*sig_x*bmd::NSIGMA*sig_y);
#endif

#if 1
	double r,theta;	/* this set is used in polar integral */
	double RR,RN;
	double dr,dtheta;
	double c,s;

	/* calculate integral using polar coordinates */
	dr = (sig_x>sig_y ? bmd::NSIGMA*sig_x/((double)(bmd::NSTEP)) : bmd::NSIGMA*sig_y/((double)(bmd::NSTEP)));
	dtheta = 2.0*TMath::Pi()/((double)(bmd::NSTEP));
	r=0.0;
	sum=0.0;
	flat=0.0;
	for( i=0; i<bmd::NSTEP; ++i )
	{
		theta = 0.0;
		for( j=0; j<bmd::NSTEP; ++j )
		{
			x = r*cos(theta);
			y = r*sin(theta);
			if( (x*x)/(bmd::NSIGMA*bmd::NSIGMA*sig_x*sig_x)+(y*y)/(bmd::NSIGMA*bmd::NSIGMA*sig_x*sig_x) < 1.0 )	/* this is three-sigma check */
			{
				rc2 = (x-off_x)*(x-off_x)+(y-off_y)*(y-off_y);
				if( rc2>testR2)
				{
					c = sig_y*r*cos(theta);
					s = sig_x*r*sin(theta);
					RR = (c*c+s*s);

					RN = RR/(2.0*sig_x*sig_x*sig_y*sig_y);
					f = exp(-RN)/(2.0*TMath::Pi()*sig_x*sig_y)*r*dr*dtheta;
					sum += f;
					flat += r*dr*dtheta;
				}
			}
			theta += dtheta;
		}
		r += dr;
	}
	flat *= (1.0-SIG3)/(M_PI*bmd::NSIGMA*sig_x*bmd::NSIGMA*sig_y);
#endif

	return sum+flat;

}

double  bmd::fpeaks(Double_t *x, Double_t *par)
{
   Double_t result = par[0] + par[1]*x[0];
   for(Int_t p=0; p<bmd::NPEAKS; ++p)
	{
      Double_t norm  = par[3*p+2];
      Double_t mean  = par[3*p+3];
      Double_t sigma = par[3*p+4];
      result += norm*TMath::Gaus(x[0],mean,sigma);
   }
   return result;
}


double bmd::profilewidth( std::vector<double>& y )
{

	int ileft = -1;
	for( int i=0; i<46; ++i )
	{
		if( y[i]<y[i+1] && y[i+1]<y[i+2] )
		{
			ileft = i;
			break;
		}
	}

	int iright = -1;
	for( int i=47; i>1; --i )
	{
		if( y[i]<y[i-1] && y[i-1]<y[i-2] )
		{
			iright = i;
			break;
		}
	}

	double width;
	if( ileft==-1 || iright==-1 )
		width = -1.0;
	else
		width = ((iright-ileft)+1.0)*0.5;

//	std::cout<<"width: "<<width<<std::endl;
	return width;
}


Int_t bmd::profilesigmaROOT( TH1F* myhist, Int_t* npeaks, Double_t* par )
{
//	fprintf(stderr,"profilesigma");
	Int_t p;
//	TCanvas* ct = new TCanvas("ct","temp",400,400);

	TH1F *hclone = (TH1F*)(myhist->Clone("hclone"));
	Double_t normh = hclone->Integral();
	hclone->Scale(1.0/normh);

//	hclone->Draw();
//	ct->Update();

	TSpectrum *s = new TSpectrum(NPEAKS);
	Int_t nfound = s->Search(hclone,3.0,"nodraw",0.002);
//	Int_t nfound = s->Search(hclone,3.0,"",0.002);
	if( nfound>3 ) nfound=3;
#if 0
	for( p=0; p<nfound; ++p )
	{
		Double_t xp = s->GetPositionX()[p];
		Double_t yp = s->GetPositionY()[p];
		printf("found peak: %f\tmean: %f\n", yp,xp);
	}
#endif
	TH1 *hb = s->Background(hclone,20,"same");

//	ct->Update();

	TF1 *fline = new TF1("fline","pol1",-12.5,12.5);
//	hb->Fit("fline");
	hb->Fit("fline","qn");

	par[0] = fline->GetParameter(0);
	par[1] = fline->GetParameter(1);
//	printf("background:\t%f\t%f\n", par[0],par[1]);

	for( p=0; p<nfound; ++p )
	{
		par[3*p+2] = s->GetPositionY()[p];
		par[3*p+3] = s->GetPositionX()[p];
		par[3*p+4] = 3;
	}

	TF1 *fit = new TF1("fit",fpeaks,-12.5,12.5,2+3*nfound);
	fit->SetParameters(par);
	fit->SetNpx(48);
	hclone->Fit("fit","qn");
//	hclone->Fit("fit");
	par[0]=fit->GetParameter(0);
	par[1]=fit->GetParameter(4);
	for( p=0; p<nfound; ++p )
	{
		par[3*p+2] = fit->GetParameter(3*p+2);
		par[3*p+3] = fit->GetParameter(3*p+3);
		par[3*p+4] = fit->GetParameter(3*p+4);
	}
	
//	printf("sigma: %f\n", par[4]);

//	ct->Update();
//	ct->WaitPrimitive();
//	ct->Close();

	*npeaks = nfound;

	return 1;
}


/**
	The following is a kludge module that translates for C++ containers to ROOT containers.
**/
int bmd::profilesigma( std::vector<double>& mwire, int* npeaks, double* par )
{
	Int_t npeaksROOT;
	Double_t parROOT[11];
	//Int_t stsROOT;
	//	int sts;
	TH1F* myhistROOT = new TH1F("myhistROOT","",48,-12.5,12.5);

	for (int i=0; i<48; ++i)
	{
		double y;
		y = mwire[i];
		myhistROOT->SetBinContent(i+1,y);
	}

	//	stsROOT =
	profilesigmaROOT( myhistROOT, &npeaksROOT, parROOT );

	//	sts = stsROOT;
	*npeaks = npeaksROOT;
	for( int i=0; i<11; ++i )
	{
		par[i] = parROOT[i];
	}

	delete myhistROOT;

	return 1;
}

double bmd::findedge( double offset, double width )
{
	double edge;

	if( offset>0.0 )
		edge =   offset + width/2.0;
	else
		edge = -(offset - width/2.0);

	return edge;
}


