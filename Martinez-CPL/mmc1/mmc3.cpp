#include "somedefs.h"
#include "somefuncs.h"
#include "simul.h"
#include <stdio.h>

int main(int argc, char *argv[])
{
    char ToFile[]="Example_001.txt"; //File to store results

    GeneralData MyData; //struct containing relevant data
    DefaultData(MyData);//setting default values

    //values for our specific example
    MyData.Mechanism=ECMECHANISM;//A+e = B k1-> <-k-1 C
    MyData.Electrode=SPHERICAL;
    MyData.R0=2;
    MyData.Gamma=1.67;
    MyData.PointsDerivs=6;
    MyData.PointsCurrent=5;
    MyData.NPointsGridSpecies=13;//points in the grid
    MyData.TimeSteps=1000;
    //bulk concentrations
    MyData.Cinfi[0]=1;//species A
    MyData.Cinfi[1]=0;//species B
    MyData.Cinfi[2]=0;//species C
    //dimensionless rate constants
    MyData.KKinec[0]=2;//k1
    MyData.KKinec[1]=0;//k-1
    //allocate memory
    double *Intensity=new double[MyData.TimeSteps];//to store Current
    double *Time=new double[MyData.TimeSteps];//to store time
    double *Potential=new double[MyData.TimeSteps];//to store Potential
    //for CV
    double PotIni,PotEnd;//start and end of scan (triangular)
    PotIni=10;
    PotEnd=-10;

    double TMax;
    //for CV
    TMax=2*fabsl(PotEnd-PotIni);//dimensionless time at the end of scan
    MyData.DeltaT=TMax/MyData.TimeSteps;

    MyData.XMax=6.0*sqrtl(TMax);

    //Calculate potential values for a CV scan
    CVScan(Potential,PotIni,PotEnd,MyData.TimeSteps);

    //objet to manage simulations with our parameters
    Simul *MySimul=new Simul(MyData);

    //Kinetic Terms to be added to Matrix
    MySimul->ObtainKineticTerms();
    //setting bulk concentrations
    MySimul->BulkConCen();
    //Filling matrix A
    MySimul->FillMatrixA(PotIni);

    //Filling Independent Terms
    MySimul->FillIndTerm();

    // analyze pattern of Matrix, only once
    MySimul->AnalyzeMatrixPattern();
    //factorize Matrix
    MySimul->FactorizeMatrix();
    //Solve linear system for the first time
    MySimul->SolveLinearSistem();
    //to obtain current
    Intensity[0]=MySimul->ObtainFlux(0);//Flux of species '0' (A)
    Time[0]=MyData.DeltaT;
    int kk;

    for(kk=1;kk<MyData.TimeSteps;kk++)
    {
        MySimul->FillIndTerm();
        //change potential
        MySimul->ChangeNernstPotential(0,MyData.NPointsGridSpecies,Potential[kk]);
        MySimul->FactorizeMatrix();//only factorize
        MySimul->SolveLinearSistem();
        Intensity[kk]=MySimul->ObtainFlux(0);
        Time[kk]=MyData.DeltaT*(double)kk;

    }

    FILE *pfile;
    pfile=fopen(ToFile,"w");
    for(kk=0;kk<MyData.TimeSteps;kk++)  //storing data into a file and showing into screen
    {
        fprintf(pfile,"%lf\t%lf\t%lf\n",Time[kk],Potential[kk],Intensity[kk]);
        printf("%lf\t%lf\t%lf\n",Time[kk],Potential[kk],Intensity[kk]);

    }
    fclose(pfile);

    //Memory free
    delete [] Time;
    delete [] Intensity;
    delete [] Potential;
    delete MySimul;

}
