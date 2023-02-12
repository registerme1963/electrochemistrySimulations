#include "somedefs.h"
#include "somefuncs.h"
#include "CoefsAlphaBeta.h"

void DefaultData(GeneralData &MData)
{
    MData.DeltaT=0.01;
    int kk;
    MData.Cinfi[0]=1;
    for(kk=1;kk<MAX_SPECIES;kk++) MData.Cinfi[kk]=0;
    MData.Gamma=1.5;
    MData.hh=0.01;
    for(kk=0;kk<MAX_STEPS_C;kk++) MData.KKinec[kk]=1;
    MData.Lambda=100;
    MData.Mechanism=EMECHANISM;
    MData.NPointsGridSpecies=12;
    MData.NSpecies=2;
    MData.PointsCurrent=5;
    MData.PointsDerivs=5;
    MData.Electrode=PLANE;
    MData.FactorG=0;
    MData.R0=2.;
    MData.TimeSteps=100;
    MData.XMax=6.;

}


double Calculatehh(GeneralData &MData)
{
    return MData.XMax*(MData.Gamma-1)/(pow(MData.Gamma,MData.NPointsGridSpecies-1)-1);
}
void CVScan(double * Perturbation,double PotIni,double PotEnd, int NPoints)
    {

    int MediumPoint;

    MediumPoint=NPoints/2;
    int kk;
    for(kk=0;kk<NPoints;kk++)
            {
            if(kk<MediumPoint)
                {
                Perturbation[kk]=PotIni+(PotEnd-PotIni)*(double)(kk)/(double)MediumPoint;
                }
            else
                {
                Perturbation[kk]=PotEnd+(PotEnd-PotIni)*(double)(MediumPoint-kk)/(double)MediumPoint;
                }
             }

    }


