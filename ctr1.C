//
// Partie 2 : Histogramme et lissage (on utilise cette foi ci TTree)
//	
void ctr1()
{

// création de l'histogramme hRes demandé dans la question 1
       TH1F *hRes = new TH1F("hRes", "Masse du Boson de Higgs",100,120.,135.);

// création des quatre histogrammes demandé dans la question 2

	// Les trois premier histogrammes: hp1x, hp1y, hp1z
       TH1F *hp1x = new TH1F("hp1x", "hp1x",100,-150.,150.);
       TH1F *hp1y = new TH1F("hp1y", "hp1y",100,-150.,150.);
       TH1F *hp1z = new TH1F("hp1z", "hp1z",100,-150.,150.);

       // l'histogramme he1
       TH1F *he1 = new TH1F("he1", "he1",100,0.,150.);

// fichier qui contient le tuple
       TFile *input = new TFile("dataNtuple.root", "read");

// création de l'objet TTree
       TTree *MyTree = (TTree*)input->Get("MyTree");

// déclaration des composantes des quadriimpulsions
       float E1, P1x, P1y, P1z, E2, P2x, P2y, P2z;

// on reli les composantes des quadriimpulsion aux branches        
       MyTree->SetBranchAddress("E1",&E1);
       MyTree->SetBranchAddress("P1x",&P1x);
       MyTree->SetBranchAddress("P1y",&P1y);
       MyTree->SetBranchAddress("P1z",&P1z);

       MyTree->SetBranchAddress("E2",&E2);
       //MyTree->SetBranchAddress("P2x",&P2x);
       //MyTree->SetBranchAddress("P2y",&P2y);
       //MyTree->SetBranchAddress("P2z",&P2z);

// nombre d'entrée disponibles
       Long64_t entries = MyTree->GetEntries();
       
// boucle sur les éléments du Tuple
       for(int i=0; i<entries; i++)
	 {
	   MyTree->GetEntry(i);
	   he1->Fill(E1);
	   hp1x->Fill(P1x);
	   hp1y->Fill(P1y);
	   hp1z->Fill(P1z);
	   
// on calcule la masse du Higgs (MHiggs= E1+E2) et on place les valeurs
//aléatoirs associée dans hRes pour les tracer à la fin
	   hRes->Fill(E1+E2);
	 }

// on cré un canvas et on le divise en quatre pads
       TCanvas *cH1 = new TCanvas("cH1","multipads",900,900);
       cH1->Divide(2,2,0,0);
       
//pour rendre mon canvas plus joli
       cH1->SetTickx();
       cH1->SetTicky();

// on trace chqu'un de nos histogrammes, he1, hp1x, hp1y et hp1z	

       cH1->cd(1);
       // les légende de l'axe de l'histogram hRes
       he1->GetXaxis()->SetTitle("E1");
       he1->Draw();

       cH1->cd(2);
       // les légende de l'axe de l'histogram hRes
       hp1x->GetXaxis()->SetTitle("P1x");
       hp1x->Draw();

       cH1->cd(3);
       // les légende de l'axe de l'histogram hRes
       hp1y->GetXaxis()->SetTitle("P1y");
       hp1y->Draw();

       cH1->cd(4);
       // les légende de l'axe de l'histogram hRes
       hp1z->GetXaxis()->SetTitle("P1z");
       hp1z->Draw();


// on crée un canvas cH2 pour tracer l'histogramme hRes
       TCanvas *cH2 = new TCanvas();
       cH2->SetTickx();
       cH2->SetTicky();
       hRes->Draw();


/*********************** Partie 3: Résultats ****************************/

       
// on peint le corp de l'istograme avec une couleur pour qu'il soit plus joli
	hRes->SetFillColor(kBlue-9);
       
// on fait un Fit (lissage) a l'histogramme hRes
       TF1 *fit = new TF1("lissage", "gaus", 0., 150.);
       hRes->Fit("lissage");

// les légendes des axes de l'histogram hRes
	hRes->GetXaxis()->SetTitle("Masse du Boson de Higge");

// Pour la masse du Higgs ainsi que l'erreur associé, voire le pdf ataché 

// pour l'encadré
       	TLegend *legend = new TLegend(0.45, 0.6, 0.78, 0.8);
	// on met la taille de la bordure égale a zéro (on l'élimine crément)
	// juste pour que notre canvas soit joli et cela grace à la commande
	// suivante 
	legend->SetBorderSize(0);
	legend->AddEntry(hRes,"les donnees", "f");
	legend->AddEntry(fit,"le fit des donnees", "l");
	legend->Draw();
}
