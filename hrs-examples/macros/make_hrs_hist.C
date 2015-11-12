#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <TFile.h>
#include <TNtuple.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TRandom3.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TMath.h>
#include <TBox.h>
#include <TPolyLine.h>
#include <TLegend.h>
void make_hrs_hist() {
  gROOT->Reset();
  gStyle->SetOptStat(0);
 gStyle->SetTitleOffset(1.,"Y");
 gStyle->SetTitleOffset(.7,"X");
 gStyle->SetLabelSize(0.04,"XY");
 gStyle->SetTitleSize(0.05,"XY");
 gStyle->SetPadLeftMargin(0.17);
     TString fstart="hrs_ran_traj_15cm_20deg";
 //   TString fstart="hrs_ran_traj_15cm_20deg";
  Double_t solang_fac=0.10*0.20*1000.;
 //   TString fstart="hrs_mod20_ran_traj_15cm_20deg";
//Double_t solang_fac=0.10*0.120*1000.;
TFile *fsnake =  new TFile("rootfiles/"+fstart+".root");
TFile *hsnake = new TFile("rootfiles/"+fstart+"_hist.root","recreate");
 Bool_t fit_matrix=kFALSE;
 Bool_t reached_fp=kTRUE;
  string newcoeffsfilename="hms_newfit_hrs_mod20_ver2.dat";
 TRandom3 err; 
//
//  string oldcoeffsfilename="hrs_recon_cosy.dat";
    string oldcoeffsfilename="hrs_newfit_standard.dat";
 // string oldcoeffsfilename="hrs_newfit_mod20.dat";
  ifstream oldcoeffsfile(oldcoeffsfilename.c_str());
  ofstream newcoeffsfile(newcoeffsfilename.c_str());
  int num_recon_terms_old;

  vector<double> xptarcoeffs_old;
  vector<double> yptarcoeffs_old;
  vector<double> ytarcoeffs_old;
  vector<double> deltacoeffs_old;
  vector<int> xfpexpon_old;
  vector<int> xpfpexpon_old;
  vector<int> yfpexpon_old;
  vector<int> ypfpexpon_old;
  vector<int> xtarexpon_old;

  vector<double> xptarcoeffs_new;
  vector<double> yptarcoeffs_new;
  vector<double> ytarcoeffs_new;
  vector<double> deltacoeffs_new;
  vector<int> xfpexpon_new;
  vector<int> xpfpexpon_new;
  vector<int> yfpexpon_new;
  vector<int> ypfpexpon_new;
  vector<int> xtarexpon_new;
  vector<double> xtartrue,ytartrue,xptartrue,yptartrue,deltatrue;
  vector<double> xfptrue,yfptrue,xpfptrue,ypfptrue;
  TString currentline;

  num_recon_terms_old = 0;

  while( currentline.ReadLine(oldcoeffsfile,kFALSE) && !currentline.BeginsWith(" ----") ){
    //    cout << currentline.Data() << endl;
    //extract the coeffs and exponents from the line:

    TString sc1(currentline(1,16));
    TString sc2(currentline(17,16));
    TString sc3(currentline(33,16));
    TString sc4(currentline(49,16));
    
    xptarcoeffs_old.push_back(sc1.Atof());
    ytarcoeffs_old.push_back(sc2.Atof());
    yptarcoeffs_old.push_back(sc3.Atof());
    deltacoeffs_old.push_back(sc4.Atof());
    
    int expontemp[5];

    for(int expon=0; expon<5; expon++){
      TString stemp(currentline(66+expon,1));
      expontemp[expon] = stemp.Atoi();
    }

    xfpexpon_old.push_back(expontemp[0]);
    xpfpexpon_old.push_back(expontemp[1]);
    yfpexpon_old.push_back(expontemp[2]);
    ypfpexpon_old.push_back(expontemp[3]);
    xtarexpon_old.push_back(expontemp[4]);

    // cout << "(C_theta, C_y, C_phi, C_delta) = (" 
    //	 << xptarcoeffs_old[num_recon_terms_old] << ", "
    //	 << ytarcoeffs_old[num_recon_terms_old] << ", " 
    //	 << yptarcoeffs_old[num_recon_terms_old] << ", "
    //	 << deltacoeffs_old[num_recon_terms_old] << "), expon = "
    //	 << xfpexpon_old[num_recon_terms_old] << xpfpexpon_old[num_recon_terms_old] 
    //	 << yfpexpon_old[num_recon_terms_old] << ypfpexpon_old[num_recon_terms_old] 
    //	 << xtarexpon_old[num_recon_terms_old] << endl;
    
    num_recon_terms_old++;
  }

  cout << "num recon terms in OLD matrix = " << num_recon_terms_old << endl;
  int npar,nfit_max=1000;
  npar=  num_recon_terms_old;
  TVectorD b_ytar(npar);
  TVectorD b_yptar(npar);
  TVectorD b_xptar(npar);
  TVectorD b_delta(npar);
  TMatrixD lambda(npar,nfit_max);

  //TMatrixD A_yptar(npar,npar);
  TMatrixD Ay(npar,npar);
 
//
 int nentries; //number of entries in file
Double_t evnum,epnum,xabs,yabs,zabs,cxabs,cyabs,czabs,mom,pathl,live,xrel,yrel,zrel,cxrel,cyrel,czrel;
 //
 TH1F *hEndPlane = new TH1F("hEndPlane"," Endplane Number",26,-.5,25.5);
 TH1F *hMom = new TH1F("hMom","Momentum (GeV)",50,.5,1.5);
 TH1F *hMomp = new TH1F("hMomp","Momentum (GeV)",50,.5,1.5);
 TH1F *hMomf = new TH1F("hMomf","Momentum (GeV)",50,.5,1.5);
 TH1F *hDeltarecon = new TH1F("hDeltarecon","Delta Recon ",40,-.12,.12);
 TH1F *hDeltadiff = new TH1F("hDeltadiff","Delta Diff % ",40,-.3,.3);
 TH1F *hDeltanew = new TH1F("hDeltanew","Delta New Recon ",40,-.12,.12);
 TH1F *hDeltanewdiff = new TH1F("hDeltanewdiff","Delta New Diff % ",40,-.3,.3);
 TH2F *hDeltadiffvsxfp = new TH2F("hDeltadiffvxtar",";Delta Diff %;xfp ",40,-.1,.1,40,-80.,80.);
 TH1F *hDelta = new TH1F("hDelta","Delta ",20,-.1,.1);
 TH1F *hDeltap = new TH1F("hDeltap","Delta ",20,-.1,.1);
 TH1F *hDeltaf = new TH1F("hDeltaf","Delta ",20,-.1,.1);
 TH1F *hDeltar = new TH1F("hDeltar",";Delta ;Solid angle (msr)  ",20,-.1,.1);
 TH1F *hytar = new TH1F("hytar","ytar (cm)",70,-7.,7.);
 TH1F *hytarrecon = new TH1F("hytarrecon","ytar recon(cm)",70,-7.,7.);
 TH1F *hytarnew = new TH1F("hytarnew","ytar new(cm)",70,-7.,7.);
 TH1F *hytardiff = new TH1F("hytardiff","ytar diff(cm)",70,-10.,10.);
 TH1F *hytarnewdiff = new TH1F("hytarnewdiff","ytar new diff(cm)",70,-10.,10.);
 TH1F *hytarp = new TH1F("hytarp","ytar (cm)",70,-7.,7.);
 TH1F *hytarf = new TH1F("hytarf","ytar (cm)",70,-7.,7.);
 TH1F *hxptar = new TH1F("hxptar","xptar ",100,-.12,.12);
 TH1F *hxptarrecon = new TH1F("hxptarrecon","xptar recon",100,-.12,.12);
 TH1F *hxptardiff = new TH1F("hxptardiff","xptar diff (mr)",100,-100,100);
 TH1F *hxptarnew = new TH1F("hxptarnew","xptar new recon",100,-.12,.12);
 TH1F *hxptarnewdiff = new TH1F("hxptarnewdiff","xptar new diff (mr)",100,-100,100);
 TH1F *hxptarp = new TH1F("hxptarp","xptar ",100,-.12,.12);
 TH1F *hxptarf = new TH1F("hxptarf","xptar ",100,-.12,.12);
 TH1F *hyptar = new TH1F("hyptar","yptar ",100,-.1,.1);
 TH1F *hyptarrecon = new TH1F("hyptarrecon","yptar recon ",100,-.1,.1);
 TH1F *hyptardiff = new TH1F("hyptardiff","yptar diff (mr) ",100,-100,100);
 TH1F *hyptarnew = new TH1F("hyptarnew","yptar new recon ",100,-.1,.1);
 TH1F *hyptarnewdiff = new TH1F("hyptarnewdiff","yptar new diff (mr) ",100,-100,100);
 TH1F *hyptarp = new TH1F("hyptarp","yptar ",100,-.1,.1);
 TH1F *hyptarf = new TH1F("hyptarf","yptar ",100,-.085,.085);
 TH2F *hxyfp = new TH2F("hxyfp","xfp vs yfp (cm)",100,-60.,60.,100,-100.,100.);
 TH2F *hxpypfp = new TH2F("hxpypfp","xpfp vs ypfp (rad)",100,-.05,.05,100,-.2,.2);
 TH2F *hxpyfp = new TH2F("hxpyfp","xpfp vs yfp ",100,-40.,40.,100,-.1,.1);
 TH2F *hxpdel = new TH2F("hxpdel","xptar vs delta ",50,-.1,.1,55,-.1,.1);
 TH2F *hxpxfp = new TH2F("hxpxfp","xfp vs xpfp ",90,-.15,.15,100,-100.,100.);
 TH2F *hypyfp = new TH2F("hypyfp","yfp vs ypfp ",100,-.06,.06,100,-10.,10.);
 TH2F *hypdel = new TH2F("hypdel","yptar vs delta ",100,-.1,.1,80,-.06,.06);
 TH2F *hypytar = new TH2F("hypytar","yptar vs ytar ",100,-10.,10.,80,-.04,.04);
 TH2F *hxpyp = new TH2F("hxpyp","xptar vs yptar ",40,-.04,.04,55,-.1,.1);
 TH2F *hxptarep = new TH2F("hxptarep","xptar vs endplane ",38,-.5,37.5,50,-.055,.055);
 TH2F *hyptarep = new TH2F("hyptarep","yptar vs endplane ",38,-.5,37.5,50,-.04,.04);
 //
 //
const int NumEndPl=26;
const int NumFocalPlane=24;
 Double_t xmin[NumEndPl]= {
,-12.,-15.,-15.,-15.,-25.
,-25.,-25.,-25.,-30.,-30.
,-30.,-30.,-30.,-30.,-30.
,-30.,-30.,-30.,-30.,-30.
,-30.,-30.,-20.,-20.,-10.,-10.
};
 Double_t xmax[NumEndPl]= {
,+12.,+15.,+15.,+15.,+25.
,+25.,+25.,+25.,+30.,+30.
,+30.,+30.,+30.,+30.,+30.
,+30.,+30.,+30.,+30.,+30.
,+30.,+30.,+20.,+20.,+10.,+10.
};
 Double_t ymin[NumEndPl]= {
,-12.,-15.,-15.,-15.,-25.
,-25.,-25.,-25.,-25.,-20.
,-30.,-30.,-30.,+50.,+40.
,-650.,-550.,-100.,+0.0,-30.
,-80.,-80.,-80.,-80.,-80.,-80.
};
 Double_t ymax[NumEndPl]= {
,+12.,+15.,+15.,+15.,+25.
,+25.,+25.,+25.,+25.,+20.
,+30.,+30.,+30.,+200.,+200.
,-550.,-450.,+100.,+100.,+30.
,+80.,+80.,+80.,+80.,+80.,+80.
};
 char *EndPlname[NumEndPl]={"Target"," 80.0cm from target"," Collimnator",
			    "Q1 ent -30.7 ","Q1 exit 47.0 "," Q1 exit","Q1Q2 61.04cm","Q1Q2 67.5cm","Q1Q2 117.0cm",
			    "Q2 91.5","Q2 center " ,"Q2 56.17","Q2 85.5","Q2D -304 cm","Q2D -150.","DIP -2681.2",
"DIP 1512cm","Dip center","DQ3 115.4 ","Q3 -57.5","Q3 91.5","Q3 exit ","Q3 57.5 ",
"Q3 85 ","Focal plane","FP +390cm"};
 TH2F *hxyEndPl[NumEndPl],*hxyEndPlp[NumEndPl],*hxyEndPlp2[NumEndPl];
 TH1F *hxEndPl[NumEndPl],*hyEndPl[NumEndPl];
 TH1F *hxEndPlp[NumEndPl],*hyEndPlp[NumEndPl];
 for ( int i = 0; i < NumEndPl ; i++) {
   hxEndPl[i] = new TH1F(Form("hxEndPl%02d", i),Form("%s; X (cm) ; Counts",EndPlname[i]),60,xmin[i],xmax[i]);
   hyEndPl[i] = new TH1F(Form("hyEndPl%02d", i),Form("%s; Y (cm) ; Counts",EndPlname[i]),60,ymin[i],ymax[i]);
   hxEndPlp[i] = new TH1F(Form("hxEndPlp%02d", i),Form("%s; X (cm) ; Counts",EndPlname[i]),60,xmin[i],xmax[i]);
   hyEndPlp[i] = new TH1F(Form("hyEndPlp%02d", i),Form("%s; Y (cm) ; Counts",EndPlname[i]),60,ymin[i],ymax[i]);
   hxyEndPl[i] = new TH2F(Form("hxyEndPl%02d", i),EndPlname[i],80,xmin[i],xmax[i],80,ymin[i],ymax[i]);
   hxyEndPlp[i] = new TH2F(Form("hxyEndPlp%02d", i),Form("%s; Y (cm) ; X (cm)  ",EndPlname[i]),80,xmin[i],xmax[i],80,ymin[i],ymax[i]);
   hxyEndPlp2[i] = new TH2F(Form("hxyEndPlp2%02d", i),Form("%s; Y (cm) ;  X (cm)",EndPlname[i]),80,xmin[i],xmax[i],80,ymin[i],ymax[i]);
 }
 //
 //
TTree *tsnake = (TTree*)fsnake->Get("ntuple");
 tsnake->SetBranchAddress("evnum",&evnum);
tsnake->SetBranchAddress("epnum",&epnum);
tsnake->SetBranchAddress("xabs",&xabs);
tsnake->SetBranchAddress("yabs",&yabs);
tsnake->SetBranchAddress("zabs",&zabs);
tsnake->SetBranchAddress("cxabs",&cxabs);
tsnake->SetBranchAddress("cyabs",&cyabs);
tsnake->SetBranchAddress("czabs",&czabs);
tsnake->SetBranchAddress("mom",&mom);
tsnake->SetBranchAddress("pathl",&pathl);
tsnake->SetBranchAddress("live",&live);
tsnake->SetBranchAddress("xrel",&xrel);
tsnake->SetBranchAddress("yrel",&yrel);
tsnake->SetBranchAddress("zrel",&zrel);
tsnake->SetBranchAddress("cxrel",&cxrel);
tsnake->SetBranchAddress("cyrel",&cyrel);
tsnake->SetBranchAddress("czrel",&czrel);
//

 Double_t xptar,yptar,xpfp,ypfp,xtar,xfp,yfp;
 int ntracks;
 Double_t cmom=0.837595,delta,ytar,ztemp,thcent;
 Double_t xsave[NumEndPl],zsave[NumEndPl];
   // In Snake abs coordinate system +x is down,+y into SHMS,+z to large angle
//
 int EndPl;
 int last_plane= NumFocalPlane;
 int nfit;
      nentries = (int)tsnake->GetEntries();
           for (int ie = 0; ie < nentries; ie++)
            {
	      //   	      cout << ie << endl;
        tsnake->GetEntry(ie);
   	  delta=(mom-cmom)/cmom;
          EndPl = epnum;
          if ( EndPl > NumEndPl) {
	    EndPl =NumEndPl;
            cout << " endplane number too large" ;
          }
	if ( epnum == 0  ) {
          ntracks++;
          yptar=czabs/cyabs;
	  xptar=cxabs/cyabs;
          ytar=zabs/10.;
          hMom->Fill(mom);
          hDelta->Fill(delta);
          hytar->Fill(ytar);
	  hxptar->Fill(xptar);
	  hyptar->Fill(yptar);
	  if (last_plane != NumFocalPlane) {
          for (int ii = 0; ii < last_plane; ii++) {
          hxyEndPlp2[ii]->Fill(zsave[ii]/10.,xsave[ii]/10.);
          }
	  }
	  last_plane=0;
        }
	  xsave[EndPl]=xrel;
	  zsave[EndPl]=zrel;
	if ( live == 0  ) {
	  if (epnum <= NumFocalPlane) last_plane=epnum;
	  hxptarep->Fill(epnum,xptar);
	   hyptarep->Fill(epnum,yptar);
	  hxptarf->Fill(xptar);
	  hyptarf->Fill(yptar);
	  hMomf->Fill(mom);
          hDeltaf->Fill(delta);
          hytarf->Fill(ytar);
          hEndPlane->Fill(epnum);
          hxyEndPl[EndPl]->Fill(zrel/10.,xrel/10.);
	}
	if ( live == 1  ) {
          hxEndPl[EndPl]->Fill(xrel/10.);
          hyEndPl[EndPl]->Fill(zrel/10.);
        }
	if ( live == 1 && epnum==NumFocalPlane   ) {
	  last_plane=NumFocalPlane;
          for (int ii = 0; ii < NumEndPl; ii++) {
          hxyEndPlp[ii]->Fill(zsave[ii]/10.,xsave[ii]/10.);
          hyEndPlp[ii]->Fill(zsave[ii]/10.);
          hxEndPlp[ii]->Fill(xsave[ii]/10.);
          }
	  hxptarp->Fill(xptar);
	  hyptarp->Fill(yptar);
          hMomp->Fill(mom);
          hDeltap->Fill(delta);
          hytarp->Fill(ytar);
	  hxyfp->Fill(zrel/10.,xrel/10.);
	  yfp=zrel/10.; // focal plane in centimeters
	  xfp=xrel/10.; // focal plane in centimeters
	  thcent=20./180.*3.14159;
          ztemp=ytar*(cos(thcent)/tan(thcent+yptar)+sin(thcent));
	  xtar=xptar*ztemp*cos(thcent)/100.;
          ypfp=czrel/cyrel;
	  xpfp=cxrel/cyrel;
	  //
	  xpfp=xpfp+err.Gaus(0.0,0.001);
	  ypfp=ypfp+err.Gaus(0.0,0.001);         
	  xfp=xfp+err.Gaus(0.0,.03);
	  yfp=yfp+err.Gaus(0.0,.03);         
	  //
	  hxpypfp->Fill(ypfp,xpfp);
	  hxpxfp->Fill(xpfp,xrel/10.);
	  hxpyfp->Fill(zrel/10.,xpfp);
	  hypyfp->Fill(ypfp,zrel/10.);
	  hxpdel->Fill(delta,xptar);
	  hxpyp->Fill(yptar,xptar);
	  hypdel->Fill(delta,yptar);
	  hypytar->Fill(ytar,yptar);
          Double_t ytartemp = 0.0,yptartemp=0.0,xptartemp=0.0,deltatemp=0.0;
          Double_t etemp;
          for( int icoeffold=0; icoeffold<num_recon_terms_old; icoeffold++ ){
        	etemp= 
	  pow( xfp / 100.0, xfpexpon_old[icoeffold] ) * 
	  pow( yfp / 100.0, yfpexpon_old[icoeffold] ) * 
	  pow( xpfp, xpfpexpon_old[icoeffold] ) * 
	  pow( ypfp, ypfpexpon_old[icoeffold] ) * 
	  pow( xtar, xtarexpon_old[icoeffold] );
        	deltatemp += deltacoeffs_old[icoeffold] * etemp;
        	ytartemp += ytarcoeffs_old[icoeffold] * etemp;
	        yptartemp += yptarcoeffs_old[icoeffold] * etemp;
	         xptartemp += xptarcoeffs_old[icoeffold] *etemp; 
	      if (nfit < nfit_max) {
              lambda[icoeffold][nfit] = etemp;
	      b_xptar[icoeffold] += (xptar) * etemp;
	      b_yptar[icoeffold] += (yptar) * etemp;
	      b_ytar[icoeffold] += (ytar) /100.0 * etemp;
	      b_delta[icoeffold] += (delta) * etemp;
              }
	  } // for loop
	  if (nfit < nfit_max) {
             nfit++;
          }
  	    xfptrue.push_back( xfp );
	    yfptrue.push_back( yfp );
	    xpfptrue.push_back( xpfp );
	    ypfptrue.push_back( ypfp );
	    xtartrue.push_back( xtar );
	    xptartrue.push_back( xptar );
	    ytartrue.push_back( ytar  );
	    yptartrue.push_back( yptar  );
	    deltatrue.push_back( delta  );
	  hytarrecon->Fill(ytartemp*100.);
	  hyptarrecon->Fill(yptartemp);
	  hxptarrecon->Fill(xptartemp);
	  hDeltarecon->Fill(deltatemp);
	  hytardiff->Fill(ytartemp*100.-ytar);
	  hyptardiff->Fill(1000*(yptartemp-yptar));
	  hxptardiff->Fill(1000*(xptartemp-xptar));
	  hDeltadiff->Fill(100*(deltatemp-delta));
	  hDeltadiffvsxfp->Fill(100*(deltatemp-delta),xfp);         
	}
	//        
	    }
// end loop over tree
  //
if ( fit_matrix) {
 for(int i=0; i<npar; i++){
    for(int j=0; j<npar; j++){
      Ay[i][j] = 0.0;
    }
 }
  for( int ifit=0; ifit<nfit; ifit++){
    if( ifit % 100 == 0 ) cout << ifit << endl;
    for( int ipar=0; ipar<npar; ipar++){
      for( int jpar=0; jpar<npar; jpar++){
	Ay[ipar][jpar] += lambda[ipar][ifit] * lambda[jpar][ifit];
      }
    }
  }
  TDecompSVD Ay_svd(Ay);
  bool ok;
  ok = Ay_svd.Solve( b_ytar );
  cout << "ytar solution ok = " << ok << endl;
  //b_ytar.Print();
  ok = Ay_svd.Solve( b_yptar );
  cout << "yptar solution ok = " << ok << endl;
  //b_yptar.Print();
  ok = Ay_svd.Solve( b_xptar );
  cout << "xptar solution ok = " << ok << endl;
  //b_xptar.Print();
  ok = Ay_svd.Solve( b_delta );
  cout << "Delta solution ok = " << ok << endl;
  //b_delta.Print();
  // calculate target quantities with new fit parameter
  for( int ifit=0; ifit<nfit; ifit++){
          Double_t ytarnew = 0.0,yptarnew=0.0,xptarnew=0.0,deltanew=0.0;
     for( int ipar=0; ipar<npar; ipar++){
       etemp=lambda[ipar][ifit];
        	deltanew += b_delta[ipar] * etemp;
        	ytarnew += b_ytar[ipar] * etemp;
	        yptarnew += b_yptar[ipar] * etemp;
	         xptarnew += b_xptar[ipar] *etemp;        
    }
	  hytarnew->Fill(ytarnew*100.);
	  hyptarnew->Fill(yptarnew);
	  hxptarnew->Fill(xptarnew);
	  hDeltanew->Fill(deltanew);
	  hytarnewdiff->Fill(ytarnew*100.-ytartrue.at(ifit));
	  hyptarnewdiff->Fill(1000*(yptarnew-yptartrue.at(ifit)));
	  hxptarnewdiff->Fill(1000*(xptarnew-xptartrue.at(ifit)));
	  hDeltanewdiff->Fill(100*(deltanew-deltatrue.at(ifit)));    
  }
  // write out coeff
  char coeffstring[100];
  cout << "writing new coeffs file" << endl;
          for( int icoeffold=0; icoeffold<num_recon_terms_old; icoeffold++ ){
      newcoeffsfile << " ";
      sprintf( coeffstring, "%16.9g", b_xptar[icoeffold] );
      newcoeffsfile << coeffstring; 
      //      newcoeffsfile << " ";
      sprintf( coeffstring, "%16.9g", b_ytar[icoeffold] );
      newcoeffsfile << coeffstring;
      sprintf( coeffstring, "%16.9g", b_yptar[icoeffold] );
      //newcoeffsfile << " ";
      newcoeffsfile << coeffstring; 
      sprintf( coeffstring, "%16.9g", b_delta[icoeffold] );
      //newcoeffsfile << " ";
      newcoeffsfile << coeffstring; 
      newcoeffsfile << " ";
	newcoeffsfile << setw(1) << setprecision(1) << xfpexpon_old[icoeffold]; 
	newcoeffsfile << setw(1) << setprecision(1) << xpfpexpon_old[icoeffold]; 
	newcoeffsfile << setw(1) << setprecision(1) << yfpexpon_old[icoeffold]; 
	newcoeffsfile << setw(1) << setprecision(1) << ypfpexpon_old[icoeffold]; 
	newcoeffsfile << setw(1) << setprecision(1) << xtarexpon_old[icoeffold]; 
      newcoeffsfile << endl;

	  }
  newcoeffsfile << " ---------------------------------------------" << endl;

  newcoeffsfile.close();
  cout << "wrote new coeffs file" << endl;
  //
	   } // if fitting
  //
    TCanvas *cep = new TCanvas("cep"," Endplanes versus xp  and yp",800,800);
    cep->Divide(1,2);
    cep->cd(1);
    hxptarep->Draw("colz");
    cep->cd(2);
    hyptarep->Draw("colz");
    //
  //
    TCanvas *c7 = new TCanvas("c7","X vs Y (cm) for endplane 1-7",800,800);
    c7->Divide(2,4);
    for (int i = 1; i <8 ; i++){
      EndPl=i;
    c7->cd(i);
    hxyEndPlp[EndPl]->Draw("box");
    //    hxyEndPlp2[EndPl]->Draw("box,same");
    hxyEndPlp2[EndPl]->SetLineColor(2);
    hxyEndPl[EndPl]->Draw("colz,same");
    }
//
    TCanvas *cq1 = new TCanvas("cq1","X vs Y (cm) for 8-13",800,800);
    cq1->Divide(2,3);
    for (int i = 1; i <7 ; i++){
      EndPl=i+7;
    cq1->cd(i);
     hxyEndPlp[EndPl]->Draw("box");
    hxyEndPl[EndPl]->Draw("colz,same");
    }
//
    TCanvas *cq2 = new TCanvas("cq2","X vs Y (cm) for 14-20",800,800);
    cq2->Divide(2,3);
    for (int i = 1; i <6 ; i++){
      EndPl=i+13;
    cq2->cd(i);
    hxyEndPlp[EndPl]->Draw("box");
    hxyEndPl[EndPl]->Draw("colz,same");
    }
//
    TCanvas *cq3 = new TCanvas("cq3","X vs Y (cm) for 20-25",800,800);
    cq3->Divide(2,3);
    for (int i = 1; i <6 ; i++){
      EndPl=i+19;
    cq3->cd(i);
    hxyEndPlp[EndPl]->Draw("box");
    hxyEndPl[EndPl]->Draw("colz,same");
    }
    TCanvas *c5 = new TCanvas("c5","Solid Angle",800,800);
    c5->Divide(1,1);
    c5->cd(1);
    hDeltar->Divide(hDeltap,hDelta,solang_fac,1.,"B");
    hDeltar->Draw();
    cout << " solid angle = " << solang_fac*(hDeltap->Integral(6,14))/(hDelta->Integral(6,14)) << " msr " << endl;
//
    TCanvas *c2 = new TCanvas("c2","Compare target",800,800);
    c2->Divide(2,2);
    c2->cd(1);
    hxptar->Draw();   
    hxptarp->SetLineColor(2);
    hxptarp->Draw("same");   
    hxptarf->SetLineColor(4);
    hxptarf->Draw("same");   
    c2->cd(2);
    hyptar->Draw();   
    hyptarp->SetLineColor(2);
    hyptarp->Draw("same");   
    hyptarf->SetLineColor(4);
    hyptarf->Draw("same");   
    c2->cd(3);
    hDelta->Draw();   
    hDeltap->SetLineColor(2);
     hDeltap->Draw("same");   
    hDeltaf->SetLineColor(4);
    hDeltaf->Draw("same");   
    c2->cd(4);
   hytar->Draw();   
    hytarp->SetLineColor(2);
     hytarp->Draw("same");   
    hytarf->SetLineColor(4);
    hytarf->Draw("same");   
//
    TCanvas *c3 = new TCanvas("c3","Focal Plane",800,800);
    TPad *c3_1 = new TPad("c3_1", "c3_1",0.01,0.51,0.49,0.99);
    TPad *c3_3 = new TPad("c3_3", "c3_3",0.51,0.51,0.99,0.99);
    TPad *c3_2 = new TPad("c3_2", "c3_2",0.01,0.01,0.49,0.49);
    TPad *c3_4 = new TPad("c3_4", "c3_4",0.51,0.01,0.99,0.49);
    c3_1->Draw();
    c3_2->Draw();
    c3_3->Draw();
    c3_4->Draw();
    c3_1->cd();
    c3_1->SetGridx(1);
    c3_1->SetGridy(1);
    hxyfp->Draw("colz");
    c3_2->cd();
    c3_2->SetGridx(1);
    c3_2->SetGridy(1);
    hxpypfp->Draw("colz");
    c3_3->cd();
    c3_3->SetGridx(1);
    c3_3->SetGridy(1);
    hxpxfp->Draw("colz");
    c3_4->cd();
    c3_4->SetGridx(1);
    c3_4->SetGridy(1);
    hypyfp->Draw("colz");
//
    TCanvas *c4 = new TCanvas("c4","2d target plots",800,800);
    c4->Divide(2,2);
    c4->cd(1);
    hxpdel->Draw("colz");
    c4->cd(2);
    hypdel->Draw("colz");
    c4->cd(3);
    hxpyp->Draw("colz");
    c4->cd(4);
    hypytar->Draw("colz");
//
    TCanvas *crecon = new TCanvas("crecon","Recon target",800,800);
    crecon->Divide(2,2);
    crecon->cd(1);
    hytarrecon->Draw();
    hytarnew->Draw("same");
    crecon->cd(2);
    hyptarrecon->Draw();
    hyptarnew->Draw("same");
    crecon->cd(3);
    hxptarrecon->Draw();
    hxptarnew->Draw("same");
    crecon->cd(4);
    hDeltarecon->Draw();
    hDeltanew->Draw("same");
    //
//
    TCanvas *cdiff = new TCanvas("cdiff","Diff target",800,800);
    cdiff->Divide(2,2);
    cdiff->cd(1);
    hytardiff->Draw();
 hytardiff->Fit("gaus");
 TF1 *fitcydiff=hytardiff->GetFunction("gaus");
    cdiff->cd(2);
    hyptardiff->Draw();
 hyptardiff->Fit("gaus");
 TF1 *fitcypdiff=hyptardiff->GetFunction("gaus");
    hyptarnewdiff->Draw("same");
    cdiff->cd(3);
    hxptardiff->Draw();
    hxptarnewdiff->Draw("same");
 hxptardiff->Fit("gaus");
 TF1 *fitcxpdiff=hxptardiff->GetFunction("gaus");
    cdiff->cd(4);
    //    hDeltadiffvsxfp->Draw("colz");
    hDeltadiff->Draw();
 hDeltadiff->Fit("gaus");
 TF1 *fitcdeldiff=hDeltadiff->GetFunction("gaus");
    hDeltanewdiff->Draw("same");
 cout << " sig_y  "<< fitcydiff->GetParameter(2) << " sig_yp " << fitcypdiff->GetParameter(2) << " sig_xp "<< fitcxpdiff->GetParameter(2) << " sig_delta " << fitcdeldiff->GetParameter(2) << endl;
    //
  hDelta->Write();
  hDeltap->Write();
  hytarp->Write();
  hyptarp->Write();
  hxptarp->Write();
  hDeltadiff->Write();
  hytardiff->Write();
  hyptardiff->Write();
  hxptardiff->Write();
 for ( int i = 0; i < NumEndPl ; i++) {
  hxEndPl[i]->Write();
  hyEndPl[i]->Write();
  hxyEndPl[i]->Write();
  hxyEndPlp[i]->Write();
  hxEndPlp[i]->Write();
  hyEndPlp[i]->Write();
  cout <<  " EndPl num = " << i << " " << EndPlname[i] << "               " << hxyEndPl[i]->Integral() << endl;
 }
   }
