void elastic(){



  //booleans
  bool just_fit=0;
  bool fresco_input = 0;
  
  //Define variables
  double D, eps_o, Z_1, Z_2, A_2, e, pi, T, diff_cs, diff_cs_mb_msr, angle_no;
  int start, end, start_temp, end_temp;
  //variables to change for different peaks
  // cout  << "25Mg elastic peaks" << endl; Z_2=12, start=1575, end=1650, A_2=25; //25Mg

  cout  << "Au elastic peaks" << endl; Z_2=79, start=1750, end=1850, A_2=197; //197Au
  //Calculate expected Rutherford cross-section
  eps_o=8.854*pow(10,-12); Z_1=1; e=1.6*pow(10,-19); pi=3.14159;
  T=8*pow(10,6)*1.6*pow(10,-19);
  D=(Z_1*Z_2*e*e)/(4*pi*eps_o*T);

  diff_cs=(D*D)/(16*pow(sin((pi*30.0/2) /180.0),4));
  // diff_cs=6.48*pow(10,-29);
  diff_cs_mb_msr=diff_cs*pow(10,31);
  diff_cs_mb_msr=diff_cs_mb_msr*pow(10,-3);

  if(fresco_input){
    // diff_cs_mb_msr=0.7070; diff_cs=0.7070*pow(10,-28); //Mg25
    diff_cs_mb_msr=28.749; diff_cs=28.749*pow(10,-28); //Au197
  }
  // diff_cs=diff_cs*pow(10,28); //converting from m^2/sr to mb/sr
  cout << "constant D is " << D << " with a diff cs in m^2/sr of " << diff_cs << " and a diff cs in mb/msr of " << diff_cs_mb_msr << endl;
  const int total_runs=107;
  const int total_angles=57;
  double run_no_array[total_runs], BCI_array[total_runs], state_counts_array[total_runs], state_value_array[total_runs], corresp_angle_array[total_runs], state_value_uncert_array[total_runs];
  double angle_array[total_angles], total_BCI_array[total_angles],total_state_value_array[total_runs],total_state_value_uncert_array[total_angles];
  int max_peak_height;
  int angle = 39;
  int run_no = 10;
  int columns = 2;
  double scaling_constant=100000;
  double deadtime, solid_angle, targ_surf_dens, integral, gaus_a,gaus_b,gaus_c,gaus_a_error,gaus_b_error, gaus_c_error, integral_error, integral_a_error, integral_c_error;
  int runs_per_angle_array[total_angles];
  int xbin,xmin,xmax,ybin,ymin,ymax;
  int peak_value; //, angle_no;
  xbin=4096;ybin=100;xmin=0;xmax=4095;ymin=0;ymax=10;


  //Get elastic histograms
  TFile *f = new TFile("all_spectra.root");

  for (run_no = 0; run_no < 107; run_no++){ //start of run loop
    // TFile *f = new TFile(Form("%ddegrees.root",angle)); //root file data is from
    if(run_no==24 || run_no==92 || run_no==93 || run_no==106){
      cout << "elastic peak at run no " << run_no << endl;
      TH1F *h1 = new TH1F("h1","h1 title", xbin,xmin,xmax);
      cout << "root input file is from run " << run_no << endl;
      if(!(f->GetListOfKeys()->Contains(Form("run%d", run_no)))){ //checking if histogram exists and if not starting again
        cout << "histogram doesn't exist: number " << run_no << endl;
        delete h1;
        continue;
      }

      //Getting histogram
      h1 = (TH1F*)f->Get(Form("run%d", run_no));
      int i, file_flag=0; //resetting flag to know if location is there or not
      int run_file_counter, corresp_angle;


      // Getting parameters for run
      // Get BCI from scaler histogram
      TFile *f_scalar = new TFile(Form("root_scalers/%d_scaler.root",run_no));
      TH1F *h_scalar=(TH1F*)f_scalar->Get("hScaler");
      //extracting scalar values
      double BCI=h_scalar->GetBinContent(1); double CLOCK=h_scalar->GetBinContent(2); double CLOCK_BUSY=h_scalar->GetBinContent(3); double MG=h_scalar->GetBinContent(10); double MG_LIVE=h_scalar->GetBinContent(11);
      deadtime=(1-CLOCK_BUSY/CLOCK);
      if(run_no<39 || run_no>63){
        solid_angle=1.0; //msr, runs 0-38, 64-106
        cout << "solid angle set to: " << solid_angle << endl;
      }
      else{
        solid_angle=1.0; //msr, runs 39-63
      }

      if (run_no>94){
        targ_surf_dens=112.0; //ug/cm^2, target no 5
      }
      else{
        targ_surf_dens=90.0; //ug/cm^2, target no 4
      }
      targ_surf_dens=(targ_surf_dens*0.000001/25)*6.022*pow(10,23); //converting to nuclei/cm^2, assuming pure 25mg (it isn't)
      //mg25:99.2, mg24:0.47, mg26:0.33
      //finding peak height
      // cout << g1->GetParameter(1) << endl;





      peak_value = h1->Integral(start,end);

      TF1 *g1 = new TF1("m1","gaus",start,end);
      h1->GetXaxis()->SetRangeUser(start,end); //zooming in
      // g1->SetParLimits(2,0,10);
      h1->Fit(g1);
      // h1->Fit(g1); //fitting
      // peak_value=g1->Integral(start,end);

      cout <<" bins counts for " << run_no << " is " << peak_value << endl;
      h1->Draw();
      g1->Draw("SAME");
      gaus_a = g1->GetParameter(0); gaus_b = g1->GetParameter(1); gaus_c = g1->GetParameter(2);


      integral=1.414*gaus_a*gaus_c*1.7725;
      cout << "integral of gaussian for " << run_no << " is " << integral << endl;
      gaus_a_error = (g1->GetParError(0)); gaus_b_error = (g1->GetParError(1)); gaus_c_error = g1->GetParError(2);
      cout << "error of parameter 0 is " << gaus_a_error  << "\nerror of parameter 1 is " << gaus_b_error << "\nerror of parameter 2 is " << gaus_c_error << endl;

      integral_a_error=gaus_c_error*1.414*1.7725*gaus_a_error;
      integral_c_error=gaus_a_error*1.414*1.7725*gaus_c_error;
      integral_error=sqrt(integral_a_error*integral_a_error + integral_c_error*integral_c_error);
      cout << "error of integral of gaussian is " << integral_error << endl;
      state_value_uncert_array[run_no]=integral_error;
      gStyle->SetOptFit();


      if(run_no==92 || run_no==106){
        cout <<"run appears to have double gaussian" << endl;
        TF1 *g1 = new TF1("m1","gaus(0)+gaus(3)",start,end);
        // TF1 *g1 = new TF1("m1","gaus",start,end);

        h1->GetXaxis()->SetRangeUser(start,end); //zooming in, needs to be before fit
        cout << "start end " << start << end << endl;
        g1->SetParameter(1,1778);          g1->SetParameter(0,9000); g1->SetParameter(2,7);
        g1->SetParameter(4,1795);        g1->SetParameter(3,9000);        g1->SetParameter(5,1);
        // g1->SetParLimits(1,1580s,1605);
        // g1->SetParLimits(2,0,8);        g1->SetParLimits(0,1,30);

        h1->Fit(g1,"0");
        // peak_value=g1->Integral(start,end);
        cout <<" bins counts for " << run_no << " is " << peak_value << endl;
        h1->Draw();
        g1->Draw("SAME");
        gaus_a = g1->GetParameter(0); gaus_b = g1->GetParameter(1); gaus_c = g1->GetParameter(2);


        integral=1.414*gaus_a*gaus_c*1.7725;
        cout << "integral of gaussian for " << run_no << " is " << integral << endl;
        gaus_a_error = (g1->GetParError(0)); gaus_b_error = (g1->GetParError(1)); gaus_c_error = g1->GetParError(2);
        cout << "error of parameter 0 is " << gaus_a_error  << "\nerror of parameter 1 is " << gaus_b_error << "\nerror of parameter 2 is " << gaus_c_error << endl;

        integral_a_error=gaus_c_error*1.414*1.7725*gaus_a_error;
        integral_c_error=gaus_a_error*1.414*1.7725*gaus_c_error;
        integral_error=sqrt(integral_a_error*integral_a_error + integral_c_error*integral_c_error);
        cout << "error of integral of gaussian is " << integral_error << endl;
        state_value_uncert_array[run_no]=integral_error;
        gStyle->SetOptFit();

      }






      //Adding values to arrays
      run_no_array[run_no]=run_no;
      BCI_array[run_no]=BCI;
      state_counts_array[run_no]=peak_value;
      if(just_fit){
        state_counts_array[run_no]=integral;
      }
      corresp_angle_array[run_no]=corresp_angle;

      cout << "BCI for run " << run_no << " is: " << BCI << endl;
      cout << "Deadtime for run " << run_no << " is: " << deadtime << endl;
      cout << "Solid angle in msr for run " << run_no << " is: " << solid_angle << endl;
      cout << "Target surface density in nuclei/cm^2 for run " << run_no << " is: " << targ_surf_dens << endl;
      cout << "Counts for run " << run_no << " is: " << peak_value << endl;
      cout << "Corresponding angle for run " << run_no << " is: " << corresp_angle  << endl;
      double total_deuts, abs_cs;
      total_deuts=BCI*(1/(1.6*pow(10,-19)))*pow(10,-10);
      cout << "number of deuterons in runs is " << total_deuts << endl;

      state_value_array[run_no]=(state_counts_array[run_no])/(total_deuts*targ_surf_dens*solid_angle*(1-deadtime));
      state_value_array[run_no]=(state_counts_array[run_no])/(total_deuts*targ_surf_dens*solid_angle*(1-deadtime));
      cout << "Cross-section value for run" << run_no << " in cm^2/msr is: " << state_value_array[run_no] << endl;
      state_value_array[run_no]=state_value_array[run_no]*pow(10,3); //converting msr to sr
      cout << "Cross-section value for run" << run_no << " in cm^2/sr is: " << state_value_array[run_no] << endl;
      state_value_array[run_no]=state_value_array[run_no]*pow(10,27);    // state_value_array[run_no]=scaling_constant*(state_value_array[run_no])/(BCI_array[run_no]);
      cout << "Cross-section value for run" << run_no << " in mb/sr is: " << state_value_array[run_no] << "" << endl;
      state_value_uncert_array[run_no]=sqrt(state_counts_array[run_no]);
      state_value_uncert_array[run_no]=state_value_uncert_array[run_no]/(total_deuts*targ_surf_dens*solid_angle*(1-deadtime));
      state_value_uncert_array[run_no]=state_value_uncert_array[run_no]*pow(10,30);
      cout << "Cross-section uncertainty value for run" << run_no << " in mb/sr is: " << state_value_uncert_array[run_no] << "" << endl;
      abs_cs=4*pi*state_value_array[run_no];
      cout << "Absolute cross-section value for run" << run_no << " in mb is: " << abs_cs << "\n" << endl;

      //Calculating surface density from known cross-section
      solid_angle=solid_angle*0.001;
      cout << "solid angle used in sr is " << solid_angle << endl;
      // targ_surf_dens=diff_cs*BCI*solid_angle*(1-deadtime)/(state_counts_array[run_no]*1.6*pow(10,-19));
      total_deuts=total_deuts/(1-deadtime);
      cout <<"deuts: " << total_deuts << " counts: " << state_counts_array[run_no]<< endl;
      // targ_surf_dens=diff_cs*total_deuts*solid_angle/(state_counts_array[run_no]*1.6*pow(10,-19));
      targ_surf_dens=(state_counts_array[run_no]*1.6*pow(10,-19))/(diff_cs*BCI*pow(10,-10)*solid_angle*(1-deadtime)); //gives target density in nuclei/m^2
      targ_surf_dens=targ_surf_dens*pow(10,-4); //converting to nuclei/cm^2
      targ_surf_dens=targ_surf_dens*A_2/(6.022*pow(10,23)); //converting to g/cm^2
      targ_surf_dens=targ_surf_dens*pow(10,6); //converting to ug/cm^2
      cout << "using known diff CS of "  << diff_cs << " in m^2/sr and a diff cs in mb/msr of " << diff_cs_mb_msr << endl;
      cout << "target density in ug/cm^2 is " << targ_surf_dens << "\n\n\n\n" << endl;
    }








    // delete h1;
  } //run loop ends

  //Calculate nuclei number density


}
