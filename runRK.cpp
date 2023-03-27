#include "functions.h"

using namespace std;

bool flagInitDone = false;


TF1 *funcPtTsallis;


void init()
{
    if( TXT_OUTPUT_TO_FILE )
    {
        // ofstream euler_output {"euler.txt"};
//        ofstream runge_kutta4_output {"runge_kutta4.txt"};
        outTxtFile = new ofstream( "runge_kutta4.txt" );
    }

    //        return 0;

    if( DO_ANIMATION )
    {
//        cout << "TEST" << endl;
        for ( int i = 0; i < num_steps+1; i++ )
            yAllSteps[i] = new StateVec( nVar );
//        cout << "TEST 2" << endl;
    }




    // #### Reid potential - https://en.wikipedia.org/wiki/Nuclear_force
    TCanvas *fCanvPotAndForce;
    if( DO_ANIMATION )
    {
        fCanvPotAndForce = new TCanvas("fCanvPotAndForce","PotAndForce",600,70,700,700);
        gPad->SetGrid();
    }
    fPotReidCoulumbAttract = new TF1( "fPotReidCoulumbAttract"
                                      , "[2]*(-10.463*TMath::Exp(-[0]*x)/[0]/x -1650.6*TMath::Exp(-4*[0]*x)/[0]/x +6484.2*TMath::Exp(-7*[0]*x)/[0]/x) - [1]/x", minPotRange, maxPotRange );

    fPotReidCoulumbAttract->SetNpx(N_POT_POINTS);
    double kC = 1.44; // k Coulomb 1.44 MeV/e^2*fm - https://en.wikipedia.org/wiki/Coulomb_constant
    double coeffStrongForce = 1;
    double strongMu = 0.7; // fm^-1, Reid
    //    fPotReidCoulumbAttract->SetParameters( 0.7, 1.44 ); // [0]=mu (fm^-1), [1] - Coulumb
    fPotReidCoulumbAttract->SetParameters( strongMu, kC, coeffStrongForce ); // [0]=mu (fm^-1), [1] - Coulumb, [2] - strong force coeff


    // ### now REPULSIVE:
    //    fPotReidCoulumbRepulsive = (TF1*)fPotReidCoulumbAttract->Clone( "fPotReidCoulumbRepulsive" );
    fPotReidCoulumbRepulsive = new TF1( "fPotReidCoulumbRepulsive"
                                        , "[2]*(-10.463*TMath::Exp(-[0]*x)/[0]/x -1650.6*TMath::Exp(-4*[0]*x)/[0]/x +6484.2*TMath::Exp(-7*[0]*x)/[0]/x) + [1]/x", minPotRange, maxPotRange );

    fPotReidCoulumbRepulsive->SetNpx(N_POT_POINTS);
    //    fPotReidCoulumbAttract->SetParameters( 0.7, 1.44 ); // [0]=mu (fm^-1), [1] - Coulumb, [2] - strong force coeff
    fPotReidCoulumbRepulsive->SetParameters( strongMu, kC, coeffStrongForce ); // [0]=mu (fm^-1), [1] - Coulumb



    // ### Reid with NO Coulomb:
    fPotReidNoCoulumb = new TF1( "fPotReidNoCoulumb", "[1]*(-10.463*TMath::Exp(-[0]*x)/[0]/x -1650.6*TMath::Exp(-4*[0]*x)/[0]/x +6484.2*TMath::Exp(-7*[0]*x)/[0]/x)", minPotRange, maxPotRange );
    fPotReidNoCoulumb->SetNpx(N_POT_POINTS);
    fPotReidNoCoulumb->SetParameters( strongMu, /*kC,*/ coeffStrongForce );

    // ### Reid with NO Coulomb:
    fPotCoulumbOnly = new TF1( "fPotCoulumbOnly", "[0]/x", minPotRange, maxPotRange );
    fPotCoulumbOnly->SetNpx(N_POT_POINTS);
    fPotCoulumbOnly->SetParameter( 0, kC );


    // ##### Av18
    TFile *filePotentialAv18 = 0x0;
//    filePotentialAv18 = new TFile( "/Users/macbookpro/cernbox/projectsIA/2020-10_Bound_state_search/toy_interactions_using_runge_kutta/PlotAv18_from_Dimitar_pp.root" );
    filePotentialAv18 = new TFile( "/Users/macbookpro/cernbox/projectsIA/2022-12-PSI_seminar_strongly_intensive/toy_interactions_using_runge_kutta/PlotAv18_from_Dimitar_pp.root" );
    TGraph *grAv18 = (TGraph*)filePotentialAv18->Get( "gAV18_1S0" );
    double nPointsAv18 = grAv18->GetN();
    double lastX_Av18 = grAv18->GetPointX(nPointsAv18-1);
    double step_x_Av18 = grAv18->GetPointX(nPointsAv18-1) - grAv18->GetPointX(nPointsAv18-2);
    int nPointsToAdd = (40 /*fm*/ -lastX_Av18) / step_x_Av18;


    TF1 *fExp = new TF1( "funcExp", "[0]*exp([1]*x)", 0, 10);
    fExp->SetParameters(-5, -2);
    fExp->SetNpx(1000);
    grAv18->Fit( "funcExp", "Q", "",  grAv18->GetPointX(nPointsAv18-40) , grAv18->GetPointX(nPointsAv18-1) );


    cout << "grAv18->GetN() = " << nPointsAv18 << endl;
    cout << "step_x_Av18 = " << step_x_Av18 << endl;
    // additional points from Reid!
    for ( int i = 0; i < nPointsToAdd; i++ )
    {
        double x = grAv18->GetPointX(nPointsAv18-1) + step_x_Av18*(i+1);
        grAv18->SetPoint( grAv18->GetN(), x, fExp->Eval(x) );
//        cout << "x = " << x << ", exp tail: " << fExp->Eval(x) << endl;

    }

    // NOW ADD COULUMB:
    TGraph *grAv18_REPULSIVE = new TGraph;
    TGraph *grAv18_ATTRACTIVE = new TGraph;
    for ( int i = 0; i < grAv18->GetN(); i++ )
    {
        double x = grAv18->GetPointX(i);
        double y = grAv18->GetPointY(i);
//        cout << "x = " << x << ", y: " << y << endl;
        grAv18_REPULSIVE->SetPoint( i, x, y + fPotCoulumbOnly->Eval(x) );
        grAv18_ATTRACTIVE->SetPoint( i, x, y - fPotCoulumbOnly->Eval(x) );

    }






    spline5_Av18 = new TSpline5( "mySplineAv18", grAv18 );
    spline5_Av18_REPULSIVE = new TSpline5( "Av18_REPULSIVE", grAv18_REPULSIVE );
    spline5_Av18_ATTRACTIVE = new TSpline5( "Av18_ATTRACTIVE", grAv18_ATTRACTIVE );


    // Now Forces:
    TGraph *gr_FORCE_Av18_REPULSIVE = new TGraph;
    TGraph *gr_FORCE_Av18_ATTRACTIVE = new TGraph;
    for ( int i = 0; i < grAv18->GetN(); i++ )
    {
        double x = grAv18->GetPointX(i);
        gr_FORCE_Av18_REPULSIVE->SetPoint( i, x, spline5_Av18_REPULSIVE->Derivative(x) );
        gr_FORCE_Av18_ATTRACTIVE->SetPoint( i, x, spline5_Av18_ATTRACTIVE->Derivative(x) );

    }


    spline5_FORCE_Av18_REPULSIVE = new TSpline5( "Av18_FORCE_REPULSIVE", gr_FORCE_Av18_REPULSIVE );
    spline5_FORCE_Av18_ATTRACTIVE = new TSpline5( "Av18_FORCE_ATTRACTIVE", gr_FORCE_Av18_ATTRACTIVE );


    //    grAv18->


    //
//    TF1 *fPotReidTEST = new TF1( "fPotReidTEST"
//                                        , "[2]*(-10.463*TMath::Exp(-[0]*x)/[0]/x -1650.6*TMath::Exp(-4*[0]*x)/[0]/x +5570*TMath::Exp(-7*[0]*x)/[0]/x) - [1]/x", minPotRange, maxPotRange );
//    fPotReidTEST->SetParameters( strongMu, kC, coeffStrongForce ); // [0]=mu (fm^-1), [1] - Coulumb

    if( DO_ANIMATION )
    {
        fPotReidCoulumbAttract->SetLineColor( kMagenta+1 );
        fPotReidCoulumbAttract->DrawCopy();

        fPotReidCoulumbRepulsive->SetLineColor( kGreen+1 );
        fPotReidCoulumbRepulsive->DrawCopy("same");



        grAv18->SetLineColor( kGray+1 );
        grAv18->SetLineWidth(3);
        grAv18->Draw("same");


        spline5_Av18->SetLineColor( kAzure );
        spline5_Av18->SetLineWidth(3);
        spline5_Av18->SetLineStyle(2);
        spline5_Av18->Draw("same");


        spline5_Av18_ATTRACTIVE->SetLineColor( kMagenta+2 );
        spline5_Av18_ATTRACTIVE->SetLineWidth(3);
        spline5_Av18_ATTRACTIVE->SetLineStyle(7);
        spline5_Av18_ATTRACTIVE->Draw("same");


        spline5_Av18_REPULSIVE->SetLineColor( kGreen+2 );
        spline5_Av18_REPULSIVE->SetLineWidth(3);
        spline5_Av18_REPULSIVE->SetLineStyle(7);
        spline5_Av18_REPULSIVE->Draw("same");




        // #### Force as derivate:
        grF_attr = new TGraphErrors;
        grF_rep = new TGraphErrors;


        for ( int i = 0; i < N_POT_POINTS; i++ )
        {
            double x = minPotRange + i*(maxPotRange - minPotRange)/N_POT_POINTS;
            grF_attr->SetPoint(i, x, fPotReidCoulumbAttract->Derivative(x) );
            grF_rep->SetPoint(i, x, fPotReidCoulumbRepulsive->Derivative(x) );
        }
        grF_attr->SetLineColor(kBlue);
//        grF_attr->Draw("same");

        grF_rep->SetLineColor(kRed+1);
//        grF_rep->Draw("same");


        // av18
//        grF_Av18 = new TGraphErrors;
//        grF_Av18_ATTRACTIVE = new TGraphErrors;
//        grF_Av18_REPULSIVE = new TGraphErrors;
//        for ( int i = 0; i < grAv18->GetN(); i++ )
//        {
//            double x = minPotRange + i*(maxPotRange - minPotRange)/N_POT_POINTS;
//            grF_Av18->SetPoint(i, x, spline5_Av18->Derivative(x) );
//            grF_Av18_ATTRACTIVE->SetPoint(i, x, spline5_Av18_ATTRACTIVE->Derivative(x) );
//            grF_Av18_REPULSIVE->SetPoint(i, x, spline5_Av18_REPULSIVE->Derivative(x) );
//        }


//        grF_Av18->SetLineColor(kGray+2);
//        grF_Av18->SetLineWidth(2);
//        grF_Av18->SetLineStyle(7);
//        grF_Av18->Draw("same");


//        grF_Av18_ATTRACTIVE->SetLineColor(kBlue);
//        grF_Av18_ATTRACTIVE->Draw("same");

//        grF_Av18_REPULSIVE->SetLineColor(kRed+1);
//        grF_Av18_REPULSIVE->Draw("same");


        fExp->DrawCopy("same");


        //
        TCanvas *fCanvQA_POTENTIAL_Av18 = new TCanvas("CanvQA_POTENTIAL_Av18","CanvQA_POTENTIAL_Av18",10,10,800,800);
        grAv18_ATTRACTIVE->SetLineColor(kBlue);
        grAv18_ATTRACTIVE->Draw("AL");
        grAv18_REPULSIVE->SetLineColor(kRed);
        grAv18_REPULSIVE->Draw("L same");
        gPad->SetGrid();


        TCanvas *fCanvQA_FORCEL_Av18 = new TCanvas("CanvQA_FORCE_Av18","CanvQA_FORCE_Av18",10,10,800,800);
        gr_FORCE_Av18_ATTRACTIVE->SetLineColor(kBlue);
        gr_FORCE_Av18_ATTRACTIVE->Draw("AL");
        gr_FORCE_Av18_REPULSIVE->SetLineColor(kRed);
        gr_FORCE_Av18_REPULSIVE->Draw("L same");
        gPad->SetGrid();


        spline5_FORCE_Av18_REPULSIVE->SetLineColor(kBlue);
        spline5_FORCE_Av18_REPULSIVE->SetLineWidth(2);
        spline5_FORCE_Av18_REPULSIVE->SetLineStyle(7);
        spline5_FORCE_Av18_REPULSIVE->Draw("same");

        spline5_FORCE_Av18_ATTRACTIVE->SetLineColor(kRed);
        spline5_FORCE_Av18_ATTRACTIVE->SetLineWidth(2);
        spline5_FORCE_Av18_ATTRACTIVE->SetLineStyle(7);
        spline5_FORCE_Av18_ATTRACTIVE->Draw("same");


    }
    //    return 0;

}

// #####
void runRKstandalone()
{
    DO_ANIMATION = 1;//false;//true;

    StateVec dummyState;
    runRK(dummyState);
}

// ############
int runRK( StateVec &outsideState ) //, bool standaloneUse = true, int nPart = -1 )
{
//    cout << "outsideState size: " << outsideState.size() << endl;

    bool standaloneUse = false;



    nPart = 20;
    if( outsideState.size() != 0 )
        nPart = outsideState.size()/6; // 6 vars per particle (3 coord, 3 mom)
    else
        standaloneUse = true;


    nVar = 6*nPart;
    StateVec initS( nVar );


//    cout << "nVar size: " << nVar << endl;

    if( standaloneUse ) // if not set from outside - set everything HERE!
    {
        funcPtTsallis = new TF1( "funcPtTsallis", "[0]*x* ([0]-1)*([0]-2)/([0]*[1])/([0]*[1]+sqrt(x*x+[2]*[2])*([0]-2)) *TMath::Power( 1+(sqrt(x*x+[2]*[2])-[2])/([0]*[1]), -[0])"
       //                                  , 0, 5 );
                                     , 0, 10 );
        // !!! for pp 7 TeV:
        funcPtTsallis->SetParameters( 11.6798, 0.247651, -0.605571 );


        for ( int iP = 0; iP < nPart; iP++ )
        {
            masses[iP] = 938;//0.2;//1;//gRandom->Uniform(0.5, 2);
            charges[iP] = 1;//gRandom->Uniform() > 0.5 ? +1 : -1;

            // pos
            if(0)
            {
                initS[6*iP+0] = iP==0 ? -2 : +2;//gRandom->Uniform(-0.7, 0.7);
                initS[6*iP+1] = iP==0 ? +0.5 : -0.5;//gRandom->Uniform(-0.7, 0.7);
                initS[6*iP+2] = 0;
            }


            double posPhi = gRandom->Uniform( 0, TMath::TwoPi() );
            if(1)
            {
                double posX, posY;
                gRandom->Rannor( posX, posY );
//                initS[6*iP+0] = 1.0*posX;//gRandom->Uniform(-2, 2);
//                initS[6*iP+1] = 1.0*posY;//gRandom->Uniform(-2, 2);
                initS[6*iP+0] = 0.5*posX;//gRandom->Uniform(-2, 2);
                initS[6*iP+1] = 0.5*posY;//gRandom->Uniform(-2, 2);
//                initS[6*iP+0] = 1.2*cos(posPhi); //1.2*posX;//gRandom->Uniform(-2, 2);
//                initS[6*iP+1] = 1.2*sin(posPhi); //1.2*posY;//gRandom->Uniform(-2, 2);
//                initS[6*iP+2] = gRandom->Uniform(-1, 1); //gRandom->Uniform(0.7, 1.7);

            }

            // mom
            if(0)
            {
                initS[6*iP+3] = iP==0 ?  400 : -400;//iP==0 ? +3 : -3;//5*gRandom->Uniform(-1, 1);
                initS[6*iP+4] = 0;
                initS[6*iP+5] = 0;
            }

            if(1)
            {
                double phiMom = posPhi;//gRandom->Uniform( posPhi - TMath::PiOver2(), posPhi + TMath::PiOver2() );
                if(phiMom<0)
                    phiMom+=TMath::TwoPi();
                else if (phiMom>TMath::TwoPi())
                    phiMom-= TMath::TwoPi();


                double pt = funcPtTsallis->GetRandom();
                initS[6*iP+3] = pt*cos(phiMom) * 1000; // 500*gRandom->Uniform(-1, 1);
                initS[6*iP+4] = pt*sin(phiMom) * 1000; // 500*gRandom->Uniform(-1, 1);
                initS[6*iP+5] = 500*gRandom->Uniform(-1, 1);//2*gRandom->Uniform(-1, 1);

            }
        }

    }
    else // WE DO HAVE SETTING OF PARTICLES FROM OUTSIDE
    {
        initS = outsideState;

        // set masses and charges
        for ( int iP = 0; iP < nPart; iP++ )
        {
            masses[iP] = 938;//0.2;//1;//gRandom->Uniform(0.5, 2);
            charges[iP] = 1;//gRandom->Uniform() > 0.5 ? +1 : -1;
        }
    }

//    StateVec outsideState(6*3);

    if( !flagInitDone )
        init();
    flagInitDone = true;


//    return 0;




//    initS[6*0+0] = 0.193479;
//    initS[6*0+1] = -2.04318;
//    initS[6*0+2] = 0.914954;
//    initS[6*0+3] = 298.023;
//    initS[6*0+4] = -102.662;
//    initS[6*0+5] = -20.5708;

//    initS[6*1+0] = -0.746627;
//    initS[6*1+1] = -0.482845;
//    initS[6*1+2] = -0.368725;
//    initS[6*1+3] = -66.2521;
//    initS[6*1+4] = -1048.11;
//    initS[6*1+5] = 12.6164;

//    initS[6*2+0] = 0.916151;
//    initS[6*2+1] = -1.14081;
//    initS[6*2+2] = -0.556665;
//    initS[6*2+3] = 583.95;
//    initS[6*2+4] = 1039.99;
//    initS[6*2+5] = 7.03905;


    if(0)for ( int iP = 0; iP < nPart; iP++ )
    {
        cout << initS[6*iP+0] << " " << initS[6*iP+1] << " " << initS[6*iP+2] << " "
                                 << initS[6*iP+3+0] << " " << initS[6*iP+3+1] << " " << initS[6*iP+3+2] << endl;

    }


    // test 2pc a-la Kepler
    if(0)
    {
        masses[0] = 10000;
        masses[1] = 1;

        charges[0] = +1;
        charges[1] = -1;

        initS[0] = 0.;
        initS[1] = 0.;
        initS[2] = 0.;
        initS[3] = 0;
        initS[4] = 0.;
        initS[5] = 0.;


        initS[6+0] = 0.;
        initS[6+1] = 0.2;
        initS[6+2] = 0.;
        initS[6+3] = -2.;
        initS[6+4] = 0.;
        initS[6+5] = 0.;
    }








    // integrate(euler_integrator, Kuramoto, print_coordinates,
    //           euler_output, initS, num_steps, dt);

    double wtime  = get_wall_time();
    double cptime = get_cpu_time();



    integrate(runge_kutta4_integrator, /*Kuramoto*/Coulomb, print_coordinates,
              /*runge_kutta4_output*/ outTxtFile, initS, num_steps, dt);


    outsideState = yFinalState;

    wtime = get_wall_time() - wtime;
    cptime = get_cpu_time() - cptime;
    if( standaloneUse )
        display_timing(wtime, cptime);



    // drawing initial pos
    if( DO_ANIMATION )
    {
        TCanvas *fCanvAnim = new TCanvas("canv_anim","Anim View",10,10,800,800);

        animRunner();

        for ( int iP = 0; iP < nPart; iP++ )
        {
            el_IS_Part[iP] = new TEllipse( 0.5 + initS[6*iP+0] / systemRadius, 0.5 + initS[6*iP+1] / systemRadius
                    , rVisNucleon, rVisNucleon );

            el_IS_Part[iP]->SetLineColor( kGray+1 );
            el_IS_Part[iP]->Draw();
        }

    }



    return 0;
}
