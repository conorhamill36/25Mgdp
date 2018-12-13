//Program to find Eex=6.26 MeV 0+ state in 26Mg, and calculate cross-sections
#include <TFile.h>
#include <TCanvas.h>
#include <TH1F.h>
// #include < >


const char* ConvertDoubleToString(double value){
    std::stringstream ss ;
    ss << value;
    const char* str = ss.str().c_str();
    return str;
}

Double_t myfunction(Double_t *x, Double_t *par){
  Float_t xx=x[0];
  Double_t f = par[0]*exp(-0.5*((xx-par[1])/par[2])*(xx-par[1])/par[2]) + par[3]*exp(-0.5*((xx-par[4])/par[5])*((xx-par[4])/par[5])); //sum of two gaussians
  return f;

}



void state_finder(){
  using namespace std;

  //booleans
  bool bg_subtract = 0;
  bool plot_angles = 1;
  bool plot_runs = 0;
  bool plot_errors = 1;
  bool set_log = 1;
  bool just_count = 1;
  bool just_fit = 0;
  bool print_runs_details = 1;
  bool print_file_details = 0;
  bool print_angle_run_correspondance = 0;

  //Declare variables
  int max_peak_height;
  int angle = 39;
  int run_no = 10;
  int columns = 2;
  int entries=4096;
  double scaling_constant=100000;
  double bg_value, bg_start, bg_end, deadtime, solid_angle, targ_surf_dens, temp_double, integral, gaus_a,gaus_b,gaus_c,gaus_a_error,gaus_b_error, gaus_c_error, integral_error, integral_a_error, integral_c_error, include_value;
  const char* temp_character;
  char parameter[2];
  const int total_runs=106;
  const int total_angles=57;
  //Setting up arrays
  double run_no_array[total_runs], BCI_array[total_runs], state_counts_array[total_runs], state_value_array[total_runs], corresp_angle_array[total_runs], state_value_uncert_array[total_runs], peak_width_array[total_runs], peak_width_uncert_array[total_runs], include_array[total_runs];
  double angle_array[total_angles], total_BCI_array[total_angles],total_state_value_array[total_runs],total_state_value_uncert_array[total_angles];
  int runs_per_angle_array[total_angles];
  int xbin,xmin,xmax,ybin,ymin,ymax;
  int peak_value, angle_no, i;
  xbin=4096;ybin=100;xmin=0;xmax=4095;ymin=0;ymax=10;
  string s;

  TFile *f = new TFile("all_spectra.root");
  TCanvas *c1 = new TCanvas();
  c1->cd();
  double starting=104;


  for (run_no = 0; run_no < total_runs; run_no++){ //start of run loop
    //Initialising arrays
    include_array[run_no]=0.0;
    run_no_array[run_no]=run_no;
    // cout << run_no_array[run_no] << " " << include_array[run_no] << endl;
  }
  ifstream peaklocationfile;
  peaklocationfile.open("2+peaklocations.txt");

  for (run_no = 0; run_no < total_runs; run_no++){ //start of run loop
    // TFile *f = new TFile(Form("%ddegrees.root",angle)); //root file data is from
    TH1F *h1 = new TH1F("h1","h1 title", xbin,xmin,xmax);
    // if(run_no==14){
    // }
    if(print_file_details){
      cout << "root input file is from run " << run_no << endl;
      }
    if(!(f->GetListOfKeys()->Contains(Form("run%d", run_no)))){ //checking if histogram exists and if not starting again
      if(print_file_details){
        cout << "histogram doesn't exist: number " << run_no << endl;
      }
      delete h1;
      continue;
    }
    //Getting histogram
    h1 = (TH1F*)f->Get(Form("run%d", run_no));
    int start, end, run_file_counter, corresp_angle;
    int file_flag=0; //resetting flag to know if location is there or not
    //looking into file to find location of peaks



    //going through input file
    for(i = 0; i < total_runs; i++){
      getline(peaklocationfile,s); //gets line and saves as string s
      const char *s_char = s.c_str(); //conversion so can be used in function
      corresp_angle=0;
      sscanf(s_char,"%d\t%d\t%d\t%d\t%lf",&run_file_counter,&start,&end, &corresp_angle, &include_value); //picking out variables
      // cout << "run no " << run_no << endl;
      cout << run_file_counter << start << end <<corresp_angle <<include_value << endl;
      if (run_file_counter == run_no){
        cout << "location of " << start << " " << end << " found for run " << run_no << endl;
        file_flag=1;
        if(include_value==1.0){
          cout << "run is being included" << endl;
          include_array[run_no]=include_value;
          // cout << include_value << endl;
          // cout << include_array[run_no] << endl;
        }
        if(include_value==0.0){
          cout << "run isn't being included" << endl;
        }
        break; // has found location data, so stops searching file
      }
    }

    if(bg_subtract){ //background subtraction
      cout <<"Background being subtracted from run " << run_no << endl;
      TH1F *hbg = new TH1F("hbg","hbg title",xbin,xmin,xmax);
      if(corresp_angle > 4 && corresp_angle < 30){
        bg_start=start-140; bg_end=start-90;
      }
      if(corresp_angle==30 || corresp_angle==33 || corresp_angle==36 || corresp_angle==39) {
        bg_start=start-200; bg_end=start-150;
      }
      if(corresp_angle > 44){
        bg_start=start+250; bg_end=end+380;
      }
      bg_value=h1->Integral(bg_start,bg_end);
      cout << "background interval is " << bg_start << " to " << bg_end << " giving a value of " <<bg_value << endl;
      bg_value=bg_value/(bg_end-bg_start);
      cout << "counts per background bin are " << bg_value << endl;

      //actual subtraction of background
      for (Int_t i=0;i<entries;i++){ //Actually filling in bins in histogram
       hbg->SetBinContent(i,bg_value);
      }
      h1->Add(hbg,-1); //taking background histogram away from used histogram
     // h1->Draw();
    }

    if (file_flag==0){ //location not found in file
      cout << "no location found for run " << run_no << " going from default positions " << endl;
      start = 1659;
      end = 1712;
    }

    // Getting parameters for run
    // Get BCI from scaler histogram
    TFile *f_scalar = new TFile(Form("root_scalers/%d_scaler.root",run_no));
    TH1F *h_scalar=(TH1F*)f_scalar->Get("hScaler");
    //extracting scalar values
    double BCI=h_scalar->GetBinContent(1); double CLOCK=h_scalar->GetBinContent(2); double CLOCK_BUSY=h_scalar->GetBinContent(3); double MG=h_scalar->GetBinContent(10); double MG_LIVE=h_scalar->GetBinContent(11);
    deadtime=(1-CLOCK_BUSY/CLOCK);
    if(run_no<39 || run_no>63){ //fixing solid angle
      solid_angle=1.0; //msr, runs 0-38, 64-106
    }
    else{
      solid_angle=0.5; //msr, runs 39-63
    }

    if (run_no>94){ //fixing surface densities
      targ_surf_dens=112.0; //ug/cm^2, target no 5
    }
    else {
      targ_surf_dens=90.0; //ug/cm^2, target no 4
    }
    targ_surf_dens=(targ_surf_dens*0.000001/25)*6.022*pow(10,23); //converting to nuclei/cm^2, assuming pure 25mg (it isn't)
    //mg25:99.2, mg24:0.47, mg26:0.33


    if (just_fit){
      if((14 < corresp_angle) && (corresp_angle<31) ){ // for obscured peaks
        cout << "run no " << run_no << " is a merged peak" << endl;
        // TF1 *g1 = new TF1("m1","gaus(0)+gaus(3)",start,1670);
        TF1 *g1 = new TF1("m1","gaus",start,end);

        h1->GetXaxis()->SetRangeUser(start,end); //zooming in, needs to be before fit
        g1->SetParameter(1,(start+end/2));          g1->SetParameter(0,1); g1->SetParameter(2,7);
        // g1->SetParameter(4,1640);        g1->SetParameter(3,1);        g1->SetParameter(5,1);
        // g1->SetParLimits(1,1580,1605);
        // g1->SetParLimits(2,0,8);        g1->SetParLimits(0,1,30);

        h1->Fit(g1,"0");
        peak_value=g1->Integral(start,end);
        peak_value=2*peak_value; //because half of peak is obscured
        cout <<" bins counts for " << run_no << " is " << peak_value << endl;
        h1->Draw();
        g1->Draw("SAME");
        gaus_a = g1->GetParameter(0); gaus_b = g1->GetParameter(1); gaus_c = g1->GetParameter(2);
        cout << "width of peak is " << gaus_c << endl;
        peak_width_array[run_no]=gaus_c;

        integral=1.414*gaus_a*gaus_c*1.7725;
        cout << "integral of gaussian for " << run_no << " is " << integral << endl;
        gaus_a_error = (g1->GetParError(0)); gaus_b_error = (g1->GetParError(1)); gaus_c_error = g1->GetParError(2);
        cout << "error of parameter 0 is " << gaus_a_error  << "\nerror of parameter 1 is " << gaus_b_error << "\nerror of parameter 2 is " << gaus_c_error << endl;
        peak_width_uncert_array[run_no] = gaus_c_error;
        integral_a_error=gaus_c_error*1.414*1.7725*gaus_a_error;
        integral_c_error=gaus_a_error*1.414*1.7725*gaus_c_error;
        integral_error=sqrt(integral_a_error*integral_a_error + integral_c_error*integral_c_error);
        cout << "error of integral of gaussian is " << integral_error << endl;
        state_value_uncert_array[run_no]=sqrt(integral);
        gStyle->SetOptFit();
      }

      else{ //non-obscured peaks
        TF1 *g1 = new TF1("m1","gaus",start,end);
        h1->GetXaxis()->SetRangeUser(start,end); //zooming in
        // g1->SetParLimits(2,0,10);
        h1->Fit(g1,"N");
        // h1->Fit(g1); //fitting
        peak_value=g1->Integral(start,end);

        cout <<" bins counts for " << run_no << " is " << peak_value << endl;
        h1->Draw();
        g1->Draw("SAME");
        gaus_a = g1->GetParameter(0); gaus_b = g1->GetParameter(1); gaus_c = g1->GetParameter(2);
        cout << "width of peak is " << gaus_c << endl;
        peak_width_array[run_no]=gaus_c;

        integral=1.414*gaus_a*gaus_c*1.7725;
        cout << "integral of gaussian for " << run_no << " is " << integral << endl;
        gaus_a_error = (g1->GetParError(0)); gaus_b_error = (g1->GetParError(1)); gaus_c_error = g1->GetParError(2);
        cout << "error of parameter 0 is " << gaus_a_error  << "\nerror of parameter 1 is " << gaus_b_error << "\nerror of parameter 2 is " << gaus_c_error << endl;
        peak_width_uncert_array[run_no] = gaus_c_error;

        integral_a_error=gaus_c_error*1.414*1.7725*gaus_a_error;
        integral_c_error=gaus_a_error*1.414*1.7725*gaus_c_error;
        integral_error=sqrt(integral_a_error*integral_a_error + integral_c_error*integral_c_error);
        cout << "error of integral of gaussian is " << integral_error << endl;
        state_value_uncert_array[run_no]=integral_error;
        gStyle->SetOptFit(); //controls what is being shown on the plot
      }
    }

    if(just_count){ //simply summing up bin content
      peak_value = h1->Integral(start,end);
      h1->GetXaxis()->SetRangeUser(start,end); //zooming in, needs to be before fit

      cout << peak_value << endl;
      state_value_uncert_array[run_no]=sqrt(peak_value);
      h1->Draw();


      if(0){ //turn on for 0+ peak
        if((14 < corresp_angle) && (corresp_angle<31) ){ // for obscured peaks

          cout << "for 0+ peak, merged peaks are being fitted" << endl;
          TF1 *g1 = new TF1("m1","gaus",start,end);

          h1->GetXaxis()->SetRangeUser(start,end); //zooming in, needs to be before fit
          // g1->SetParameter(1,(start+end/2));          g1->SetParameter(0,1);
          // g1->SetParameter(2,7);
          // g1->SetParameter(4,1640);        g1->SetParameter(3,1);        g1->SetParameter(5,1);
          // g1->SetParLimits(1,1580,1605);
          g1->SetParameter(2,5.5);
          g1->SetParLimits(2,4,8);    //    g1->SetParLimits(0,1,30);
          h1->Fit(g1,"0");
          peak_value=g1->Integral(start,end);
          cout <<" bins counts for " << run_no << " is " << peak_value << endl;
          if (run_no==95){
            h1->Draw();
            g1->Draw("SAME");
          }

          gaus_a = g1->GetParameter(0); gaus_b = g1->GetParameter(1); gaus_c = g1->GetParameter(2);
          cout << "width of peak is " << gaus_c << endl;
          peak_width_array[run_no]=gaus_c;

          integral=1.414*gaus_a*gaus_c*1.7725;
          peak_value=integral;
          cout << "integral of gaussian for " << run_no << " is " << integral << endl;
          gaus_a_error = (g1->GetParError(0)); gaus_b_error = (g1->GetParError(1)); gaus_c_error = g1->GetParError(2);
          cout << "error of parameter 0 is " << gaus_a_error  << "\nerror of parameter 1 is " << gaus_b_error << "\nerror of parameter 2 is " << gaus_c_error << endl;
          peak_width_uncert_array[run_no] = gaus_c_error;
          integral_a_error=gaus_c_error*1.414*1.7725*gaus_a_error;
          integral_c_error=gaus_a_error*1.414*1.7725*gaus_c_error;
          integral_error=sqrt(integral_a_error*integral_a_error + integral_c_error*integral_c_error);
          cout << "error of integral of gaussian is " << integral_error << endl;
          state_value_uncert_array[run_no]=sqrt(integral);
          gStyle->SetOptFit();
        }
      }
    }


    if (run_no==24 || run_no==86 || run_no==92 || run_no==93 || run_no==106){ //excluding elastic runs
      cout << "run_no " << run_no << " is an elastic run so counts set to zero " << endl;
      state_value_uncert_array[run_no]=0;
      peak_value=0;
      corresp_angle=0;
    }


    //Adding values to arrays
    run_no_array[run_no]=run_no;
    BCI_array[run_no]=BCI;
    state_counts_array[run_no]=peak_value;
    corresp_angle_array[run_no]=corresp_angle;
    double total_deuts;
    total_deuts=BCI*(1.0/(1.6*pow(10,-19)))*pow(10,-10);

    if(print_runs_details){
      cout << "BCI for run " << run_no << " is: " << BCI << endl;
      cout << "Deadtime for run " << run_no << " is: " << deadtime << endl;
      cout << "Solid angle in msr for run " << run_no << " is: " << solid_angle << endl;
      cout << "Target surface density in nuclei/cm^2 for run " << run_no << " is: " << targ_surf_dens << endl;
      cout << "Counts for run " << run_no << " is: " << peak_value << endl;
      cout << "Corresponding angle for run " << run_no << " is: " << corresp_angle  << endl;
      cout << "number of deuterons in runs is " << total_deuts << endl;
      cout << "is file being included? " << include_value << endl;
    }

    state_value_array[run_no]=(state_counts_array[run_no])/(total_deuts*targ_surf_dens*solid_angle*(1-deadtime));
    // state_value_array[run_no]=(state_counts_array[run_no])/(total_deuts*targ_surf_dens*solid_angle*(1-deadtime));
    // cout << "Cross-section value for run" << run_no << " in cm^2/msr is: " << state_value_array[run_no] << endl;
    state_value_array[run_no]=state_value_array[run_no]*pow(10,3);
    // cout << "Cross-section value for run" << run_no << " in cm^2/sr is: " << state_value_array[run_no] << endl;
    state_value_array[run_no]=state_value_array[run_no]*pow(10,27);    // state_value_array[run_no]=scaling_constant*(state_value_array[run_no])/(BCI_array[run_no]);
    cout << "Cross-section value for run" << run_no << " in mb/sr is: " << state_value_array[run_no] << endl;
    state_value_uncert_array[run_no]=state_value_uncert_array[run_no]/(total_deuts*targ_surf_dens*solid_angle*(1-deadtime));
    state_value_uncert_array[run_no]=state_value_uncert_array[run_no]*pow(10,30);
    cout << "Cross-section uncertainty value for run" << run_no << " in mb/sr is: " << state_value_uncert_array[run_no] << "\n" << endl;


    // delete h1;
  } //run loop ends
  peaklocationfile.close();


  cout <<"Printing include array \n\n\n\n\n\n\n" << endl;
  for (i=0;i < total_runs; i++){
    cout << run_no_array[i] << " include?" << include_array[i] << endl;
  }


  //angles being looped through
  int angle_flag=0;
  int runs_per_angle=0;
  cout << "\n\n\nRun loop ended" << endl << "angle loop starts\n\n" << endl;
  for (angle_no = 0; angle_no < total_angles; angle_no++){ //angle loops begins
    // cout << "\nangle number is " << angle_no << endl;
    angle_array[angle_no]=angle_no;
    angle_flag=0;
    runs_per_angle=0;
    // cout << "run nested loop begins for angle " << angle_no << endl;
    for (run_no = 0; run_no < total_runs; run_no++){ //run nested loop begins
      if(corresp_angle_array[run_no]==angle_no){ //matching value found
        // cout << "run number " << run_no << " matches angle " << angle_no << endl;
        // cout << "with a cross-sectin of " << state_value_array[run_no] << "and an uncertainty of " << state_value_uncert_array[run_no] << endl;
        // cout << include_array[run_no] << endl;
          // cout << "BCI for run " << run_no << " is " << BCI_array[run_no] << "" << endl;

        if(include_array[run_no]==1.0){
          // cout << "input file says to include run no " <<  run_no << " for angle number " << angle_no << endl;
          total_BCI_array[angle_no]=total_BCI_array[angle_no]+BCI_array[run_no];
          //adding values of integrals together
          total_state_value_array[angle_no]+=state_value_array[run_no];
          //adding uncertainties together
          total_state_value_uncert_array[angle_no]+=((state_value_uncert_array[run_no])*(state_value_uncert_array[run_no]));
          angle_flag=1;
          runs_per_angle++;
        }

        if(include_array[run_no]==0.0){
          // cout << "input file says not to include run no " <<  run_no << " for angle number " << angle_no << endl;
        }

      }
    }//run nested loop ends

    if(angle_flag==0){
      total_state_value_array[angle_no]=0;
    }

    runs_per_angle_array[angle_no]=runs_per_angle;
    if(runs_per_angle>0){
      total_state_value_array[angle_no]=(total_state_value_array[angle_no])/(runs_per_angle);
      // total_state_value_uncert_array[angle_no]=(sqrt(state_value_uncert_array[angle_no]));
      starting=total_state_value_uncert_array[angle_no];
      total_state_value_uncert_array[angle_no]=(sqrt(starting))/(runs_per_angle);
      // total_state_value_uncert_array[angle_no]=state_value_uncert_array[angle_no]/(runs_per_angle);
      cout << "for angle " << angle_no << " value is " << total_state_value_array[angle_no] << " uncert is " << total_state_value_uncert_array[angle_no] << " with " << runs_per_angle << " runs per angle "<< endl;
    }
    // cout << "total BCI for angle " << angle_no << " is " << total_BCI_array[angle_no] << endl;

  } //angle loop ends


  for (angle_no = 1; angle_no < total_angles; angle_no++){
    // total_state_value_uncert_array[angle_no]=total_state_value_uncert_array[angle_no];
    if (total_state_value_array[angle_no] > 0.00000000000001){
      // cout << "arrays: " << angle_array[angle_no] << " " << total_state_value_array[angle_no] << " " << total_state_value_uncert_array[angle_no] << " " << " " << total_BCI_array[angle_no] << " " << runs_per_angle_array[angle_no] << endl;
    }
    if (total_state_value_array[angle_no] < 0.00000000000001 || total_state_value_array[angle_no] > 1000000000000000){
      total_state_value_array[angle_no]=0.0;
    }

  }


  for (run_no = 0; run_no < total_runs; run_no++){ // resetting arrays because of ROOT's short term memory loss
    run_no_array[run_no]=run_no;

    if(state_value_array[run_no] < 0.0000000001 || state_value_array[run_no] > 100000000000.0){
      state_value_array[run_no] = 0.0;
    }

    if(state_value_uncert_array[run_no] < 0.0000000001 || state_value_uncert_array[run_no] > 100000000000.0 || run_no==80){
      state_value_uncert_array[run_no] = 0.0;
    }

    if(run_no > 39 && run_no < 91){
      state_value_uncert_array[run_no] = 0.0;
      state_value_array[run_no] = 0.0;
    }

    if(peak_width_array[run_no] < 0.0000000000001 || peak_width_array[run_no] > 1000000000000000){
      peak_width_array[run_no]=0.0;
    }

    if(print_runs_details){
      cout << "arrays: " << run_no_array[run_no] << " " <<state_value_array[run_no] << " " << state_value_uncert_array[run_no] << endl;
      // cout << "arrays: " << run_no_array[run_no] << " " << peak_width_array[run_no] << endl;
    }

    if(corresp_angle_array[run_no]==0){
      state_value_array[run_no] = 0.0;
      state_value_uncert_array[run_no] = 0.0;
    }
    // state_value_uncert_array[run_no]=0.1;
    // cout << "arrays: " << run_no << " " << state_value_array[run_no] << endl;
  }

  total_state_value_array[0]=0;//fudge because stupid things are happening with zeroes




  //Plotting begins

  TCanvas *c2 = new TCanvas("c2","Graph",200,10,500,300);

  if(plot_runs){
    cout << "runs being plotted " << endl;
    TGraphErrors* gr = new TGraphErrors(total_runs,run_no_array,state_value_array,0,state_value_uncert_array);

    if(plot_errors){
      // TGraphErrors* gr = new TGraphErrors(total_runs,run_no_array,state_value_array,0,state_value_uncert_array);

      cout << "errors being plotted " << endl;
    }
    else{
    }
    // TGraph* gr = new TGraph(total_runs,run_no_array,peak_width_array);
    gr->GetXaxis()->SetRange(0,total_runs);

    for(run_no=0; run_no < total_runs; run_no++){//adding angle labels to points
      if(corresp_angle_array[run_no] < 0.01){
        continue;
      }
      temp_character=ConvertDoubleToString(corresp_angle_array[run_no]);
      TLatex *latex = new TLatex(gr->GetX()[run_no], gr->GetY()[run_no],temp_character);
      latex->SetTextSize(0.02);
      gr->GetListOfFunctions()->Add(latex);
    }

    gr->GetXaxis()->SetTitle("Run no"); gr->GetYaxis()->SetTitle("#frac{d#sigma}{d#Omega} (mb/sr)");
    // gr->GetYaxis()->SetRange(0.001,1.0);
    gr->Draw("A*");
  }

  total_state_value_uncert_array[11]=0.0;

  total_state_value_uncert_array[47]=0.0;

  total_state_value_uncert_array[41]=0.0;


  if(plot_angles){
    cout << "angles being plotted " << endl;
    // if(plot_errors) {
    TGraphErrors* gr = new TGraphErrors(total_angles,angle_array,total_state_value_array,0,total_state_value_uncert_array);
      // gr->GetXaxis()->SetTitle("angle/#circ"); gr->GetYaxis()->SetTitle("#frac{d#sigma}{d#Omega} (mb/sr)");
      // gr->Draw("A*");
    // }
    // else {
      // TGraph* gr = new TGraph(total_angles,angle_array,total_state_value_array);
    gr->GetXaxis()->SetTitle("angle/#circ"); gr->GetYaxis()->SetTitle("#frac{d#sigma}{d#Omega} (mb/sr)");
    gr->Draw("A*P");
    // }
  }

  if(set_log){
    c2->SetLogy();
  }
  c2->SetGrid();

  // f->Close();
  return;
}
