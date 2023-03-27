#include "functions.h"



/*------------------------------------------------------------*/
StateVec Coulomb (const StateVec &x)
{

    StateVec f;
    f = StateVec( 6*nPart );
#pragma omp parallel for reduction(+:sumx)


    bool flag_all_dist_increased = true;
    bool flag_all_dist_above_limit = true;

    for ( int i=0; i < nPart; i++ ) // particle loop
    {
        if(PRINT_DEBUG) cout << "##### STARTING PARTICLE " << i << endl;
        // x[i + 0,1,2] - coordinates
        // x[i + 3,4,5] - velocities

        // derivatives for coordinates: velocities
        double p2 = 0;
        for ( int d = 0; d < 3; d++ )  // dimensions loop
        {
            double p_d = x[6*i + 3+d];
            p2 += p_d*p_d;
        }
        double m = masses[i];
        double E = sqrt(m*m + p2);
        for ( int d = 0; d < 3; d++ )  // dimensions loop
        {
            double p_d = x[6*i + 3+d];
            double beta = p_d/E;
            if(0 && i==0) cout << "beta = " << beta << endl;

            f[6*i+d] = beta;
        }


        // derivatives for momenta: dp/dt=F
        double Fproj[3] = {};
        for(int j=0; j < nPart; j++ )
        {
            if (i == j)
                continue;

            double dr2 = 0;
            double dx[3];
            for ( int d = 0; d < 3; d++ )  // dimensions loop
            {
                //                dx[d] = x[i+d] - x[j+d];
                dx[d] = x[6*i+d] - x[6*j+d];
                if(PRINT_DEBUG) cout << "dx[d] = " << dx[d] << endl;
                dr2 += dx[d]*dx[d];
            }
            if(PRINT_DEBUG) cout << "dr2 = " << dr2 << endl;
            //            double F = 6*0.5 * charges[i]*charges[j] /dr2;
            double r = sqrt(dr2);
            double F = 0;
            if( r>0.05 && r<40 ) // r limits by hand!
            {
                double dEdr = 0; // the -Force
//                if( charges[i]*charges[j] < 0 ) // attraction
//                    dEdr = 1*( r>minPotRange ? /*fPotReidCoulumbAttract*/spline5_Av18_ATTRACTIVE->Derivative(r) : /*fPotReidCoulumbAttract*/spline5_Av18_ATTRACTIVE->Derivative(minPotRange+0.01) );
//                else // repulsion
//                    dEdr = 1*(r>minPotRange ? /*fPotReidCoulumbRepulsive*/spline5_Av18_REPULSIVE->Derivative(r) : /*fPotReidCoulumbRepulsive*/spline5_Av18_REPULSIVE->Derivative(minPotRange+0.01) );
                double rMod = r>0.05 ? r : 0.05; // TAKE LOWER EDGE BY HAND!!!
                if( charges[i]*charges[j] < 0 ) // attraction
                    dEdr = spline5_FORCE_Av18_ATTRACTIVE->Eval(rMod);
                else // repulsion
                    dEdr = spline5_FORCE_Av18_REPULSIVE->Eval(rMod);

                F = -dEdr; // * charges[i]*charges[j]; // /dr2;
            }


            double cosTheta = dx[2]/r; // z/r
            double sinTheta = sqrt(1 - cosTheta*cosTheta);
            double rxy = r*sinTheta;

            double cosPhi = dx[0]/rxy;   // x/r
            double sinPhi = sqrt(1 - cosPhi*cosPhi);
            if(PRINT_DEBUG) cout << "cosTheta = " << cosTheta << ", sinTheta = " << sinTheta << ", cosPhi = " << cosPhi << ", sinPhi = " << sinPhi << endl;

            F *= 0.00001; // SUPPRESS BY HAND!!
            Fproj[0] += F*sinTheta*cosPhi;
            Fproj[1] += F*sinTheta*sinPhi * (dx[1]>0 ? 1 : -1);
            Fproj[2] += F*cosTheta;


            // check if distances grow
            if( global_flag_check_distances )
            {
                if( r < prev_dr[i][j] )
                    flag_all_dist_increased = false;
                prev_dr[i][j] = r;
                if( r < MAX_DIST_BETWEEN_PARTICLES )
                    flag_all_dist_above_limit = false;
            }


        }

        // now assign accelerations:
        for ( int d = 0; d < 3; d++ )  // dimensions loop
        {
            f[6*i + 3+d] = Fproj[d];// / masses[i];
            if(PRINT_DEBUG) cout << "i = " << i << ", F[" << d << "] = " << Fproj[d] << endl;
        }
    }

    // QA:
    if( 0 && global_flag_check_distances )
    {
        cout << "flag_all_dist_increased = " << flag_all_dist_increased << ", flag_all_dist_above_limit = " << flag_all_dist_above_limit << endl;
    }

    // check termination condition
    if( global_flag_check_distances && flag_all_dist_increased && flag_all_dist_above_limit )
    {
//        global_terminate_flag = true;
    }


    // QA check:
    for ( int i = 0; i < 6; i++ )  // QA
    {
        if(PRINT_DEBUG) cout << "AFTER i = " << i << ", f[" << i << "] = " << f[i] << endl;
    }

    return f;
}
/*------------------------------------------------------------*/
void print_coordinates (double time, const StateVec &y, std::ofstream &output) {
    output <<std::setprecision(10);
    //    int N = PARN;
    int nVar = 6*nPart;
    for (int i=0; i<nVar; i++)
    {
        //        output <<std::setw(20) << sin(y[i]);
        output <<std::setw(20) << y[i];
    }
    output<< std::endl;
}
/*------------------------------------------------------------*/

void euler_integrator (StateVec &y, DerivFunc dydt, double dt) {
    y += dydt (y) * dt;
}

/*------------------------------------------------------------*/
void runge_kutta4_integrator (StateVec &y, DerivFunc dydt, double dt) {
    global_flag_check_distances = true;
    auto k1 = dydt (y);
    global_flag_check_distances = false;
    auto k2 = dydt (y + dt*k1/2.);
    auto k3 = dydt (y + dt*k2/2.);
    auto k4 = dydt (y + dt*k3);
    y += (k1 + 2.*k2 + 2.*k3 + k4) * dt/6;


    // QA check:
    if(PRINT_DEBUG) for ( int i = 0; i < 6*nPart; i++ )  // QA
    {
        cout << "y[" << i << "] = " << y[i] << endl;
    }




    //visualize "nucleons"
    //    if(0)for ( int iP = 0; iP < nPart; iP++ )
    //    {
    //        elPart[iP][timeStepCounter] = new TEllipse( 0.5 + fVisNucleusRadiusNucleus * y[6*iP+0] / systemRadius, 0.5 + fVisNucleusRadiusNucleus * y[6*iP+1] / systemRadius
    //                , rVisNucleon, rVisNucleon );

    //    }


    //    timeStepCounter++;
}

// ##### my animation
int global_step = 0;
//StateVec global_y;
//double global_dt;
//OutputFunc global_output_func;
//std::ofstream global_output;
//DerivFunc global_dydt;
//Integrator global_integrator;
//int global_num_steps;


TTimer *timer;

//Float_t phiPad = 30;



/*------------------------------------------------------------*/
void integrate (Integrator integrator, DerivFunc dydt, OutputFunc output_func, 
                std::ofstream *output, StateVec y, int num_steps, double dt)
{
    for (auto step = 0; step < num_steps; ++step)
    {
        // remember initial state:
        if( DO_ANIMATION && step==0 )
            *(yAllSteps[0]) = y;

        // QA check full mom conservation:
        if(0)
        {
            double tot_px = 0;
            double tot_py = 0;
            double tot_pz = 0;
            for ( int iP = 0; iP < nPart; iP++ )
            {
                tot_px += y[6*iP+3];
                tot_py += y[6*iP+4];
                tot_pz += y[6*iP+5];
            }
            cout << "tot_px = " << tot_px << ", tot_py = " << tot_py << ", tot_pz = " << tot_pz << endl;
        }


        if( TXT_OUTPUT_TO_FILE )
            output_func (step*dt, y, *output);
        integrator (y, dydt, dt);



        if( DO_ANIMATION )
            *(yAllSteps[step+1]) = y;

        if( global_terminate_flag )
        {
            global_terminate_flag = false;
            for(int i=0; i < nPart; i++ )
                for(int j=0; j < nPart; j++ )
                {
                    prev_dr[i][j] = 0;
                }

//            cout << "terminated at step " << step << " out of " << num_steps << endl;
            break;
        }
    }
    yFinalState = y;
}
/*------------------------------------------------------------*/
double get_wall_time(){
    struct timeval time;
    if (gettimeofday(&time,NULL)){
        //  Handle error
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}
/*------------------------------------------------------------*/
double get_cpu_time(){
    return (double)clock() / CLOCKS_PER_SEC;
}
/*------------------------------------------------------------*/
void display_timing(double wtime, double cptime)
{
    using namespace std;
    int wh,ch;
    int wmin,cpmin;
    double wsec,csec;
    wh = (int) wtime/3600;
    ch = (int) cptime/3600;
    wmin = ((int)wtime % 3600)/60;
    cpmin = ((int)cptime % 3600)/60;
    wsec = wtime-(3600.*wh+60.*wmin);
    csec = cptime-(3600.*ch+60.*cpmin);
    printf ("Wall Time : %d hours and %d minutes and %.4f seconds.\n",wh,wmin,wsec);
    printf ("CPU  Time : %d hours and %d minutes and %.4f seconds.\n",ch,cpmin,csec);
}
/*------------------------------------------------------------*/






void animRunner()
{

    for ( int iP = 0; iP < nPart; iP++ )
    {
        double factorR = masses[iP] < 2. ? masses[iP]  : 2; // artificially limit the R
        elAnimation[iP] = new TEllipse( 0, 0, rVisNucleon*factorR, rVisNucleon*factorR );
        int color = charges[iP] > 0 ? kRed+1 : kBlue+1;
        elAnimation[iP]->SetLineColorAlpha( color, 0.9 );
        elAnimation[iP]->SetFillColorAlpha( color, 0.1 );

        double pmPos[3] = {};
        pmAnimation[iP] = new TPolyMarker3D( 1, pmPos, 20 );
        pmAnimation[iP]->SetMarkerColorAlpha( color, 0.8 );
        pmAnimation[iP]->SetMarkerSize( factorR );


        for ( int j = 0; j < nPart; j++ )
        {
            int lineId = iP*nPart+j;
            int lineColor = (charges[iP] == charges[j] ) ? kBlue+1 : kRed+1;
            lineAnimation[lineId] = new TLine(0,0,1,1);
            lineAnimation[lineId]->SetLineColorAlpha( color, 0.2 );

            double plPos[6] = {};
            plAnimation[lineId] = new TPolyLine3D( 2, plPos );
            plAnimation[lineId]->SetLineColorAlpha( color, 1 );//0.2 );
        }
    }



    timer = new TTimer(50);//20);
    timer->SetCommand("Animate()");
    timer->TurnOn();

}

void Animate()
{
    StateVec yThis = *yAllSteps[global_step];

    double dist12; // dist to monitor!
    double kT12; // dist to monitor!
    {
        int pId1 = 0;
        int pId2 = 1;
        double x1 = yThis[6*pId1+0];
        double y1 = yThis[6*pId1+1];
        double z1 = yThis[6*pId1+2];
        double x2 = yThis[6*pId2+0];
        double y2 = yThis[6*pId2+1];
        double z2 = yThis[6*pId2+2];
        double dx = x2-x1;
        double dy = y2-y1;
        double dz = z2-z1;
        dist12 = sqrt( dx*dx + dy*dy + dz*dz );

        double px1 = yThis[6*pId1+3];
        double py1 = yThis[6*pId1+4];
        double pz1 = yThis[6*pId1+5];
        double px2 = yThis[6*pId2+3];
        double py2 = yThis[6*pId2+4];
        double pz2 = yThis[6*pId2+5];
        double dpx = px2-px1;
        double dpy = py2-py1;
        double dpz = pz2-pz1;
        kT12 = sqrt( dpx*dpx + dpy*dpy + dpz*dpz );

    }
    cout << "doing step " << global_step << "..." << "  dist12 = " << dist12 << ", kT12 = " << kT12 << endl;
    //just in case the canvas has been deleted
    if (!gROOT->GetListOfCanvases()->FindObject("canv_anim")) return;
    //    t += step;
    //   f2->SetParameter(0,TMath::Cos(t));


    // ##### MAIN CALC Runge-Kutta
    //    global_output_func (global_step*global_dt, global_y, &global_output);
    //    global_integrator (global_y, global_dydt, global_dt);




    // ##### drawing
    //    gPad->Clear();
    if(1)for ( int iP = 0; iP < nPart; iP++ )
    {
        elAnimation[iP]->SetX1( 0.5 +  yThis[6*iP+0] / systemRadius );
        elAnimation[iP]->SetY1( 0.5 +  yThis[6*iP+1] / systemRadius );
        elAnimation[iP]->Draw();

        if(0)
        {
            double x =  yThis[6*iP+0] / systemRadius;
            double y =  yThis[6*iP+1] / systemRadius;
            double z =  yThis[6*iP+2] / systemRadius;
            pmAnimation[iP]->SetPoint(0, x, y, z );
            //            cout << pmAnimation[iP]->GetN() << endl;
            pmAnimation[iP]->Draw( "p" );
        }


        // 2D or 3D lines
        if(0)for ( int j = 0; j < nPart; j++ )
        {
            int lineId = iP*nPart+j;
            double x1 =  yThis[6*iP+0] / systemRadius;
            double y1 =  yThis[6*iP+1] / systemRadius;
            double z1 =  yThis[6*iP+2] / systemRadius;

            double x2 =  yThis[6*j+0] / systemRadius;
            double y2 =  yThis[6*j+1] / systemRadius;
            double z2 =  yThis[6*j+2] / systemRadius;

            double dx = yThis[6*j+0] - yThis[6*iP+0];
            double dy = yThis[6*j+1] - yThis[6*iP+1];
            double dz = yThis[6*j+2] - yThis[6*iP+2];
            double dist = sqrt( dx*dx + dy*dy + dz*dz );

            lineAnimation[lineId]->SetX1( 0.5 + x1 );
            lineAnimation[lineId]->SetX2( 0.5 + x2 );
            lineAnimation[lineId]->SetY1( 0.5 + y1 );
            lineAnimation[lineId]->SetY2( 0.5 + y2 );

            plAnimation[lineId]->SetPoint( 0, x1, y1, z1 );
            plAnimation[lineId]->SetPoint( 1, x2, y2, z2 );

            if( dist < 0.8 )
            {
                double width = 2./dist;
                if ( width > 8 )
                    width = 8;
                lineAnimation[lineId]->SetLineWidth( width );
                plAnimation[lineId]->SetLineWidth( width );
                //                lineAnimation[lineId]->Draw();
                plAnimation[lineId]->Draw();
            }
        }



    }


    //    phiPad += 2;
    //    gPad->SetPhi(phiPad);
    gPad->Modified();
    gPad->Update();

    global_step += 1;
    if(global_step == num_steps)
        timer->TurnOff();
}

