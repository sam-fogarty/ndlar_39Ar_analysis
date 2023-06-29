#include "TFile.h"
#include "TH1D.h"
#include "TH1.h"
#include <iostream>
#include <fstream>
using namespace std;

// program to read the radioactive decay TGraphs from
// $LARSOFT_DATA_DIR/Radionuclides on DUNE and turns them into
// histograms. Then produces an ascii text file to use in
// a geant4 macro.
void TGraphToHisto(){
  TFile *root_File = TFile::Open("/sdf/home/s/sfogarty/Desktop/RadDecay/analysis/radiological_spectra/ROOT/Potassium_40.root");
  TGraph* graph = 0;
  root_File->GetObject("Gammas", graph);

  int nbins = 100;
  cout << "nbins = " << nbins << endl;
  auto nPoints = graph->GetN(); // get # of points in graph
  cout << "Number of Points in Graph = " << nPoints << endl;

  // create variables for x and y points
  double x,y;
  double x_last, y_last, x_start, y_start;
  graph->GetPoint(nPoints-1,x_last,y_last); // final point
  graph->GetPoint(0, x_start,y_start); // starting point
  double binsize = (x_last - x_start)/(float)nbins;
  cout << "Binsize = " << binsize << " keV" << endl;

  TH1D *energies = new TH1D("energies", "Decay Energies Histogram",nbins,x_start,x_last);

  double start = 0.0;
  double end = binsize;
  double ysum = 0.0;

  double totaly = 0.0;
  for (int j=0; j<nPoints;j++){
    graph->GetPoint(j,x,y);
    totaly = totaly + y;
  }

  // loop thru points in graph, fill histogram
  for (int j=0;j<nPoints;j++)
    {
      graph->GetPoint(j,x,y);
      energies->Fill(x, y/totaly);
    }

  double mean = energies->GetMean();
  cout << "Mean energy = " << mean << " keV" << endl;

  TCanvas* energiesCanvas = new TCanvas("energiesCanvas", "energies",640,480);
  energies->Draw("hist");
  energiesCanvas->Print("K40GammaDecayEnergiesHisto.pdf","Title:40K Gamma Decay Energies");
  //cout << "Number of bins = " << h->GetNbinsX() << endl;

  auto c1 = new TCanvas("c1","Graph of Decay Energies vs Rate",700,500);
  c1->SetGrid();
  graph->SetMarkerColor(4);
  graph->SetMarkerStyle(21);
  graph->Draw();
  c1->Print("K40GammaDecayEnergiesGraph.pdf","Title:Decay Energies");
  // convert histogram to ascii file for geant4 macro
  ofstream ene_txt("Potassium_40_Gammas.txt");
  //ene_txt << "f.write('/gps/hist/point " << 0.0 << " " << 0.0 << "\\n')" << endl;
  if (ene_txt.is_open())
    {
      for (int i=0;i<nbins;i++)
	{
	  ene_txt << "f.write('/gps/hist/point " << to_string(binsize*(float)(i+1)*0.001) << " " << energies->GetBinContent(i) << "\\n')" << endl;

	}
//f.write("/gps/hist/point 0.028250 0.000000\n")
    } else {cout << "The file failed to open." << endl;}
  ene_txt.close();


}
