#include <iostream>
#include <TFile.h>
#include <TSystem.h>
#include <TTree.h>
#include <TCint.h>
#include <TNtuple.h>
void make_snake_ntuple() {
   TString dir = gSystem->WorkingDirectory();
   cout << "file = " << dir.Data() << endl;
   TFile *f = new TFile(Form("%s/rootfiles/hrs_ran_traj_15cm_20deg.root",dir.Data()),"RECREATE");
   TTree *T = new TTree("ntuple","snake data");
   Long64_t nlines = T->ReadFile(Form("%s/fort.31",dir.Data()),"evnum/D:epnum/D:xabs:yabs:zabs:cxabs:cyabs:czabs:mom:pathl:live:xrel:yrel:zrel:cxrel:cyrel:czrel");
   printf(" found %lld points\n",nlines);
   cout << "file out = " << f << endl;
   T->Write();
   //   T->Close();
}
