//
// Ce code calcule les quardiimpulsions des deux photons qui provient de
// la désintegration de la résonnance Boson de Higgs (on utilise ici TNtuple) et les stoque dans
// dataNtuple.root
//

void ntuple()
{ 
// on utilise la classe TRANDOM3
      auto gen = new TRandom3();

// on déclare les composantes des deux quadriimpulsions des deux photons
      double P1x, P1y, P1z, P2x, P2y, P2z, E1, E2;

// on crée le fichier  dataNtuple.root qui va contenir les données 
      TFile *output = new TFile("dataNtuple.root", "recreate");

      TNtuple *MyTree = new TNtuple("MyTree", "MyTree", "E1:P1x:P1y:P1z:E2:P2x:P2y:P2z");


  for(int i =0; i<100000; i++)
    {
// on génère la masse du Higgs on utilisant la loi de Breit-Wigner en fesant
// attention car d'après ce qu'on a vue en TSD cette loi peut donner des 
// valeurs négatives!! or la masse des particules n'est pas négative donc on 
// doit éliminer les valeurs négatives
      Double_t MHiggs;
       do{
      MHiggs = gen->BreitWigner(125.0,0.4);
        }while(MHiggs<0);
       
// on tenant conpte de l'erreur qui provient de la résolution du détecteur
// MHiggs devient
      Double_t sigma = 0.01*MHiggs;
      MHiggs += gen->Gaus(0.,sigma);

// on génère aléatoirement Theta et Phi selon une loi Uniforme
      Double_t Phi,Theta;
      Phi= gen->Uniform(0.,TMath::TwoPi());
      Theta = gen->Uniform(0.,TMath::Pi());
      

// on écrit les expressions des composantes des deux quadriimpulsions
      
      E1=MHiggs/2. ;
      P1x=(MHiggs/2.)*cos(Phi)*sin(Theta) ;
      P1y=(MHiggs/2.)*sin(Phi)*sin(Theta);
      P1z=(MHiggs/2.)*cos(Theta);

      E2=MHiggs/2 ;
      P2x=-(MHiggs/2)*cos(Phi)*sin(Theta) ;
      P2y=-(MHiggs/2)*sin(Phi)*sin(Theta);
      P2z=-(MHiggs/2)*cos(Theta);
      
      MyTree->Fill(E1,P1x,P1y,P1z,E2,P2x,P2y,P2z);
    }
  
// ecriture dans un fichier dataNtuple.root
      output->Write();

// on ferme notre fichier
      output->Close();
  
}
