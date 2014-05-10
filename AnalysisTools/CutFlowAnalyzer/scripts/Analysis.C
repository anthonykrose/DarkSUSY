void Analysis(){


  //    TFile *f = new TFile("darksusy_mD04_ctau05mm_test.root");
  //    TFile *f = new TFile("darksusy_mD04_ctau05mm_test2.root");
  //  TFile *f = new TFile("darksusy_mD0250.root");
  //  TFile *f = new TFile("darksusy_mD0250_ctauExp_02.root");
  //  TFile *f = new TFile("darksusy_mD0250_ctauExp_2.root");
  //  TFile *f = new TFile("darksusy_mD0250_ctauExp_5.root");
   TFile *f = new TFile("darksusy_mD0550_ctauExp_5.root");
//   TFile *f = new TFile("darksusy_mD0700_ctauExp_5.root");


  f->cd("Analysis");
  
  f->ls();


  Int_t ev;
  Int_t   genm;
  Int_t   trig;
  Float_t ptgenm[20];
  Float_t etagenm[20];

  Int_t recm;
  Float_t ptrecm[20];
  Float_t etarecm[20];

  Int_t mjets;
  Int_t mj1mu;
  Int_t mj2mu;

  Float_t ptmumj1[20];
  Float_t ptmumj2[20];

  Float_t etamumj1[20];
  Float_t etamumj2[20];

  Float_t mj1dz;
  Float_t mj2dz;

  Float_t mj1mass;
  Float_t mj2mass;

  Float_t mj1iso;
  Float_t mj2iso;

  Int_t vtx;

  TTree *t = (TTree*)f->Get("Analysis/Events");


  t->SetBranchAddress("event",&ev);
  t->SetBranchAddress("trigger",&trig);
  t->SetBranchAddress("isVtx",&vtx);

  t->SetBranchAddress("genmuons",&genm);
  t->SetBranchAddress("ptgenMuons",&ptgenm);
  t->SetBranchAddress("etagenMuons",&etagenm);
  
  t->SetBranchAddress("recmuons",&recm);
  t->SetBranchAddress("ptrecMuons",&ptrecm);
  t->SetBranchAddress("etarecMuons",&etarecm);

  t->SetBranchAddress("recmujets",&mjets);
  t->SetBranchAddress("mj1muons",&mj1mu);
  t->SetBranchAddress("mj2muons",&mj2mu);

  t->SetBranchAddress("ptmj1muons",&ptmumj1);
  t->SetBranchAddress("ptmj2muons",&ptmumj2);
  t->SetBranchAddress("etamj1muons",&etamumj1);
  t->SetBranchAddress("etamj2muons",&etamumj2);

  t->SetBranchAddress("dzmj1",&mj1dz);
  t->SetBranchAddress("dzmj2",&mj2dz);

  t->SetBranchAddress("massmj1",&mj1mass);
  t->SetBranchAddress("massmj2",&mj2mass);

  t->SetBranchAddress("isomj1",&mj1iso);
  t->SetBranchAddress("isomj2",&mj2iso);

  Int_t nentries = t->GetEntries();

  Int_t mcount[4]={0};
  Int_t rmcount[4]={0};
  Int_t mjcountb=0;
  Int_t mjcount=0;
  Int_t mjcount2mu=0;
  Int_t mjcountdz=0;
  Int_t mjcounttrig=0;
  Int_t mjcountm=0;
  Int_t mjcountiso=0;
  Int_t mjcountvtx=0;
  
  for(int k=0;k<nentries;k++){
    
    t->GetEntry(k);

    Int_t hptm_count=0;
    Int_t m_count=0;
    for(int j=0;j<genm;j++){
      if(ptgenm[j]>17.0 && fabs(etagenm[j])<0.9) hptm_count++;
      if(ptgenm[j]>8.0 && fabs(etagenm[j])<2.4) m_count++;
    }
    
    if(hptm_count>0){
      mcount[0]++;
      if(m_count>1) mcount[1]++;
      if(m_count>2) mcount[2]++;
      if(m_count>3) mcount[3]++;
    }


    Int_t hptrm_count=0;
    Int_t rm_count=0;
    for(int j=0;j<recm;j++){
      if(ptrecm[j]>17.0 && fabs(etarecm[j])<0.9) hptrm_count++;
      if(ptrecm[j]>8.0 && fabs(etarecm[j])<2.4) rm_count++;
    }
    
    if(hptrm_count>0){
      rmcount[0]++;
      if(rm_count>1) rmcount[1]++;
      if(rm_count>2) rmcount[2]++;
      if(rm_count>3) rmcount[3]++;
    }

    if(mjets==2){
     
      
      mjcountb++;

      bool hptmj1=false;
      bool hptmj2=false;

      //      cout<<" muons in mujet1   "<<mj1mu<<endl;
      //      cout<<" muons in mujet2   "<<mj2mu<<endl;

      for(int m=0;m<mj1mu;m++){
	if(ptmumj1[m]>17.0 && fabs(etamumj1[m])<0.9){ hptmj1=true;break;}
      }
      for(int m=0;m<mj2mu;m++){
	if(ptmumj2[m]>17.0 && fabs(etamumj2[m])<0.9){ hptmj2=true;break;}
      }

      //      cout<<" event  "<<ev<<"  ptmj1  "<<ptmumj1[0]<<" etamj1  "<<etamumj1[0]<<"  ptmj2  "<<ptmumj2[0]<<" etamj2  "<<etamumj2[0]<<"  bool1  "<<hptmj1<<" bool2  "<<hptmj2<<endl;


      if(hptmj1 || hptmj2){ mjcount++;

	//      cout<<" event  "<<ev<<"  ptmj1  "<<ptmumj1[0]<<" etamj1  "<<etamumj1[0]<<"  ptmj2  "<<ptmumj2[0]<<" etamj2  "<<etamumj2[0]<<"  bool1  "<<hptmj1<<" bool2  "<<hptmj2<<endl;

	double dz_cut = 0.1;
	if((hptmj1||hptmj2) && (mj1mu==2&&mj2mu==2)){
	  mjcount2mu++;
// 	  std::cout << "nominal_no_bs_VertexDzDiff->Fill(" << fabs(mj1dz - mj2dz) << ");" << std::endl;
	  if(fabs(mj1dz - mj2dz)<dz_cut) mjcountdz++;
	  if(fabs(mj1dz - mj2dz)<dz_cut && trig==1) mjcounttrig++;
	  if(fabs(mj1dz - mj2dz)<dz_cut && trig==1 &&  fabs(mj1mass-mj2mass) < (0.13 + 0.065*(mj1mass+mj2mass)/2.0)) mjcountm++;
	  if(fabs(mj1dz - mj2dz)<dz_cut && trig==1 &&  fabs(mj1mass-mj2mass) < (0.13 + 0.065*(mj1mass+mj2mass)/2.0) && 
	     (mj1iso<2.0 && mj2iso<2.0)) mjcountiso++;
	  if(fabs(mj1dz - mj2dz)<dz_cut && trig==1 &&  fabs(mj1mass-mj2mass) < (0.13 + 0.065*(mj1mass+mj2mass)/2.0) && 
	     (mj1iso<2.0 && mj2iso<2.0) && vtx==1) mjcountvtx++;
	}
      }
    }
  }
  
  cout<<"  Total Events          "<<t->GetEntries()<<"    Rel. Eff    "<<t->GetEntries()/t->GetEntries()<<endl;
  cout<<"  pT>17.0 && |eta|<0.9  "<<mcount[0]<<"    Rel. Eff       "<<mcount[0]/(t->GetEntries()*1.0)<<endl;
  cout<<"  pT>8.0 && |eta|<2.4   "<<mcount[1]<<"    Rel. Eff       "<<mcount[1]/(mcount[0]*1.0)<<endl;
  cout<<"  pT>8.0 && |eta|<2.4   "<<mcount[2]<<"    Rel. Eff       "<<mcount[2]/(mcount[1]*1.0)<<endl;
  cout<<"  pT>8.0 && |eta|<2.4   "<<mcount[3]<<"    Rel. Eff       "<<mcount[3]/(mcount[2]*1.0)<<endl;
  cout<<"  Acceptance (gen)                "<<mcount[3]/(t->GetEntries()*1.0)<<endl;
  cout<<"  pT>17.0 && |eta|<0.9  "<<rmcount[0]<<"    Rel. Eff       "<<rmcount[0]/(t->GetEntries()*1.0)<<endl;
  cout<<"  pT>8.0 && |eta|<2.4   "<<rmcount[1]<<"    Rel. Eff       "<<rmcount[1]/(rmcount[0]*1.0)<<endl;
  cout<<"  pT>8.0 && |eta|<2.4   "<<rmcount[2]<<"    Rel. Eff       "<<rmcount[2]/(rmcount[1]*1.0)<<endl;
  cout<<"  pT>8.0 && |eta|<2.4   "<<rmcount[3]<<"    Rel. Eff       "<<rmcount[3]/(rmcount[2]*1.0)<<endl;
  cout<<"  Acceptance (rec)                "<<rmcount[3]/(t->GetEntries()*1.0)<<endl;
  cout<<"  muJets==2 (bef pt,eta)  "<<mjcountb<<"    Rel. Eff       "<<mjcountb/(rmcount[3]*1.0)<<endl;
  cout<<"  muJets==2               "<<mjcount<<"    Rel. Eff       "<<mjcount/(rmcount[3]*1.0)<<endl;
  cout<<"  2mu per jet             "<<mjcount2mu<<"    Rel. Eff       "<<mjcount2mu/(mjcount*1.0)<<endl;
  cout<<"  |dz|<1mm                "<<mjcountdz<<"    Rel. Eff       "<<mjcountdz/(mjcount2mu*1.0)<<endl;
  cout<<"  trigger                 "<<mjcounttrig<<"    Rel. Eff       "<<mjcounttrig/(mjcountdz*1.0)<<endl;
  cout<<"  mass1 ~ mass2           "<<mjcountm<<"    Rel. Eff       "<<mjcountm/(mjcounttrig*1.0)<<endl;
  cout<<"  isolation  2            "<<mjcountiso<<"    Rel. Eff       "<<mjcountiso/(mjcountm*1.0)<<endl;
  cout<<"  vtx                     "<<mjcountvtx<<"    Rel. Eff       "<<mjcountvtx/(mjcountiso*1.0)<<endl;
  cout<<"  Auxiliary efficiency    "<<mjcountvtx/(rmcount[3]*1.0)<<endl;
  cout<<"  Full efficiency         "<<mjcountvtx/(t->GetEntries()*1.0)<<endl;
  cout<<"  Full efficiency/gen       "<<( (rmcount[3]/(t->GetEntries()*1.0))*mjcountiso/(mjcount*1.0)) /(mcount[3]/(t->GetEntries()*1.0)) <<endl;


  


}

//  LocalWords:  mjcountvtx
