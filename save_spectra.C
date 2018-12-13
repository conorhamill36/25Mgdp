//Program takes in all .txt files and compiles them into histograms in a single root file
//C. B. Hamill

void save_spectra(){
  cout << "starting saving spectra" << endl;
  //Declare variables
  int angle = 39;  int run_no;  int columns = 2; double total_count;
  const int total_runs = 107;
  char parameter[2]; double run_array[total_runs], peak_height_array[total_runs], total_count_array[total_runs];

  //Variables to maybe change
  gROOT->SetBatch(kTRUE);
  // TFile *f_output = new TFile("all_spectra.root","RECREATE"); //where all spectra will be saved to
  for (run_no = 0; run_no < total_runs; run_no++){ //start of loop going over runs
    // TFile *f = new TFile(Form("%ddegrees.root",angle)); //root file data is from

    //Checking if the file exists.
    const char *file = Form("run%d.txt.root",run_no);
    ifstream ifile(file);
    if (!ifile.is_open()) {
      printf("File for run %d does not exist \tStarting loop again.\n",run_no);
      run_array[run_no] = 0;
      peak_height_array[run_no] = 0;
      continue; // starts loop again
    }
    ifile.close();

    TFile *f = new TFile(Form("run%d.txt.root",run_no)); //root file data is from
    // cout << "root input file is from run " << run_no << endl;
    // declaring pointer to ntuple
    TNtuple *myNTuple;
    //new myNTuple is ntuple from file
    myNTuple = (TNtuple*)(f->Get("ntuple"))->Clone("myNTuple");
    const int entries = myNTuple->GetEntries();
    const int no_entries = entries; //gets around weird root rules

    Float_t p[no_entries], position[no_entries], counts[no_entries];

    for (Int_t i=0;i<columns;i++){
      sprintf(parameter,"p%d",i);
      myNTuple->SetBranchAddress(parameter,&p[i]); //called when cloning trees
    }

    for (Int_t i=0;i<entries;i++){
      myNTuple->GetEntry(i);
      position[i]=p[0];
      counts[i]=p[1];
      // cout << p[0] << "  " << p[1] << endl;
    }

    //making histogram
    int xbin,xmin,xmax,ybin,ymin,ymax;
    xbin=4096;ybin=100;xmin=0;xmax=4095;ymin=0;ymax=10;


    //Allows looking at spectra in histogram plots
    TCanvas *c1 = new TCanvas();
    TH1F *hPosition = new TH1F("hPosition",Form("run%d.txt.root",run_no),xbin,xmin,xmax);
    hPosition->GetXaxis()->SetTitle("Position (ch)");
    hPosition->GetYaxis()->SetTitle("Counts");

    //filling in histogram
    for (Int_t i=0;i<entries;i++){ //Actually filling in bins in histogram
      // hPosition->Fill(position[i],counts[i]);
      for(Int_t j=0; j < counts[i]; j++){
        hPosition->Fill(position[i]);
      }
      // hPosition->SetBinContent(position[i],counts[i]);
    }

    // hPosition->Draw("hist");
    // gFile = f_output;
    TFile *f_output = new TFile("all_spectra.root","UPDATE"); //where all spectra will be saved to
    hPosition->Write(Form("run%d",run_no));
    f_output->Close();

    double max_peak_height, new_height;
    int i, file_flag=0; //resetting flag to know if location is there or not
    int start, end, run_file_counter;


  } //run loop ends


  return;
}
