#include "simul.h"
#include "somefuncs.h"
#include "CoefsAlphaBeta.h"

Simul::Simul(GeneralData Data)
{
    MyData=Data;
    switch(MyData.Mechanism)
        {
        case EMECHANISM:MyData.NSpecies=2;break;
        case ECMECHANISM:MyData.NSpecies=3;break;//other cases can be added
        default: MyData.NSpecies=2;
        }
    switch(MyData.Electrode)
    {
        case PLANE:MyData.FactorG=0.0;break;
        case CYLINDRICAL:MyData.FactorG=1.0;break;
        case SPHERICAL:MyData.FactorG=2.0;break;
    }

    MyData.hh=Calculatehh(MyData);
    MyData.Lambda=MyData.DeltaT/MyData.hh/MyData.hh;

    MyData.OneRow=MyData.NPointsGridSpecies*MyData.NSpecies;
    Gammai=new double[MyData.NPointsGridSpecies];
    Gamma2i=new double[MyData.NPointsGridSpecies];
    XX=new double[MyData.NPointsGridSpecies];


    ObtainAlphaAndBeta();

    int kk;

    for(kk=0;kk<MyData.NPointsGridSpecies;kk++)
        {
        Gammai[kk]=pow(MyData.Gamma,(double)kk);
        Gamma2i[kk]=Gammai[kk]*Gammai[kk];
        }

    XX[0]=0;
    for(kk=1;kk<MyData.NPointsGridSpecies;kk++)
        {
        XX[kk]=XX[kk-1]+MyData.hh*Gammai[kk-1];
        }

    int NzPerRow,NzTotal;
    NzPerRow=MyData.PointsDerivs*2; //aproximate value (in excess) of Non Zeros per row
    NzTotal=NzPerRow*MyData.OneRow;

    Concentration=new double[MyData.OneRow];
    TermIndp=new double[MyData.OneRow];
    MContainer.reserve(NzTotal);
    MatrixA.resize(MyData.OneRow,MyData.OneRow);
    MatrixA.reserve(NzTotal);
    Vectorb.resize(MyData.OneRow);
    Vectorx.resize(MyData.OneRow);


}
Simul::~Simul()
{
    delete [] Gammai;
    delete [] Gamma2i;
    delete [] XX;
    delete [] Concentration;
    delete [] TermIndp;
}

void Simul::FillMatrixA(double Potential)
{

    //Sparse Matrix is filled using Triplets (row,column,value)
    MatrixA.setZero();
    SurfaceConditions(Potential);

    int jj,ii,mm;
    int iiMat;
    double value;
    int EndNormal;
    EndNormal=MyData.NPointsGridSpecies-(MyData.PointsDerivs-2);

    for(mm=0;mm<MyData.NSpecies;mm++)
        {

        for(ii=1;ii<MyData.NPointsGridSpecies;ii++)
            {
            if(ii<EndNormal)
               for(jj=-1;jj<MyData.PointsDerivs-1;jj++)
                    {
                    iiMat=ii+mm*MyData.NPointsGridSpecies;//row in MatrixA
                    value=CoeffMatrixN2(ii,jj);
                    MContainer.push_back(MTriplet(iiMat,iiMat+jj,value));
                    }

            else
               for(jj=-1;jj<MyData.NPointsGridSpecies-ii;jj++)
                    {
                    iiMat=ii+mm*MyData.NPointsGridSpecies;//row in MatrixA
                    value=CoeffMatrixN2(ii,jj);
                    MContainer.push_back(MTriplet(iiMat,iiMat+jj,value));

                    }

            }
        }
    AddKineticToMatrix();
    //effective filling of matrix, repeat terms are added
    MatrixA.setFromTriplets(MContainer.begin(), MContainer.end());

}
void Simul::FillIndTerm()
{
    int ii,mm;

    for(mm=0;mm<MyData.NSpecies;mm++)
    {
        //Zero at matrix rows storing surface conditions
        TermIndp[mm*MyData.NPointsGridSpecies]=0;

        for(ii=1;ii<MyData.NPointsGridSpecies;ii++)
           TermIndp[ii+mm*MyData.NPointsGridSpecies]=ObtainBN2(ii,mm);
    }


}

double Simul::ObtainBN2(int ii,int IndexSpecies)
{
    double SumGeom;
    SumGeom=MyData.FactorG*MyData.hh*Gammai[ii]/(MyData.R0+XX[ii]);
    double value;
    int EndNormal;
    EndNormal=MyData.NPointsGridSpecies-(MyData.PointsDerivs-2);
    int aux,aux2;
    aux=ii-EndNormal;
    aux2=MyData.PointsDerivs-2;

    int jj;
    double ToAdd=0;
    if(ii>=EndNormal)
            {
            for(jj=0;jj<=aux;jj++)
                    ToAdd+=(CoefBeta[aux2-jj+1]*SumGeom+CoefAlpha[aux2-jj+1])*MyData.Cinfi[IndexSpecies];

            }
    value=-Concentration[ii+MyData.NPointsGridSpecies*IndexSpecies]*Gamma2i[ii]/MyData.Lambda-ToAdd;

    return value;

}

void Simul::ObtainKineticTerms()
{


    //by default, no kinetic terms
    int kk;
    for (kk=0;kk<MAX_SPECIES;kk++)TKinO1[kk].TotalTerms=0;
    switch(MyData.Mechanism)
    {

    case EMECHANISM:
        return;//nothing to do

        //EC Mechanism:  A + e = B k1-> <-k-1 C
    case ECMECHANISM:


        TKinO1[1].TotalTerms=2;//Species 1 (B), two contributions
        TKinO1[1].Index[0]=1;//contribution of species 1 (B)
        TKinO1[1].Index[1]=2;//contribution of species 2 (C)
        TKinO1[1].Terms[0]= -MyData.KKinec[0];//k1, value multipying [B]
        TKinO1[1].Terms[1]= MyData.KKinec[1];//k-1, value multiplying [C]

        TKinO1[2].TotalTerms=2;//Species 2 (C), two contributions
        TKinO1[2].Index[0]=1;//contribution of species 1 (B)
        TKinO1[2].Index[1]=2;//contribution of species 2 (C)
        TKinO1[2].Terms[0]= MyData.KKinec[0];//k1, value multipying [B]
        TKinO1[2].Terms[1]= -MyData.KKinec[1];//k-1, value multiplying [C]

        return;
//other cases can be added

    }


}

void Simul::AddKineticToMatrix()
{
    //no problems with duplicate entries in MContainer, they are added after
    int mm,kk,ii,iiMat,jjMat;
    double hkin;
    double hh2;
    hh2=MyData.hh*MyData.hh;
    for(mm=0;mm<MyData.NSpecies;mm++)
    {
        for(kk=0;kk<TKinO1[mm].TotalTerms;kk++)
        {
            for(ii=1;ii<MyData.NPointsGridSpecies;ii++)
            {

                hkin=hh2*Gamma2i[ii];
                iiMat=ii+mm*MyData.NPointsGridSpecies;
                jjMat=ii+TKinO1[mm].Index[kk]*MyData.NPointsGridSpecies;
                MContainer.push_back(MTriplet(iiMat,jjMat,TKinO1[mm].Terms[kk]*hkin));
            }
        }
    }


}
void Simul::BulkConCen()
{
    int kk,jj;
    for(kk=0;kk<MyData.NPointsGridSpecies;kk++)
    {
        for(jj=0;jj<MyData.NSpecies;jj++)
            Concentration[jj*MyData.NPointsGridSpecies+kk]=MyData.Cinfi[jj];
    }
}

void Simul::ObtainAlphaAndBeta()
{
    int kk;
    for(kk=0;kk<MyData.PointsCurrent;kk++)CoefBeta0[kk]=Beta_N_1(MyData.PointsCurrent,kk,MyData.Gamma);

    for(kk=0;kk<MyData.PointsDerivs;kk++)
    {
        CoefAlpha[kk]=Alpha_N_2(MyData.PointsDerivs,kk-1,MyData.Gamma);
        CoefBeta[kk]=Beta_N_2(MyData.PointsDerivs,kk-1,MyData.Gamma);
    }
}
double Simul::CoeffMatrixN2(int ii, int PosRel)
{
    double SumGeom;
    SumGeom=MyData.FactorG*MyData.hh*Gammai[ii]/(MyData.R0+XX[ii]);

    switch(PosRel)
    {
        case -1:return CoefAlpha[0]+SumGeom*CoefBeta[0];
        case 0: return CoefAlpha[1]-Gamma2i[ii]/MyData.Lambda+SumGeom*CoefBeta[1];
        default: return CoefAlpha[PosRel+1]+SumGeom*CoefBeta[PosRel+1];
    }


}
void Simul::SolveLinearSistem()
{
    int kk;
    for(kk=0;kk<MyData.OneRow;kk++)Vectorb[kk]=TermIndp[kk];
    Vectorx=MSolver.solve(Vectorb);
    for(kk=0;kk<MyData.OneRow;kk++)Concentration[kk]=Vectorx[kk];


}
double Simul::ObtainFlux(int Species)
{

        int kk;
        double Add;
        Add=0;
        for(kk=0;kk<MyData.PointsCurrent;kk++)
            {
            Add+=CoefBeta0[kk]*Concentration[Species*MyData.NPointsGridSpecies+kk];
            }
        return Add/MyData.hh;

}
void Simul::ChangeNernstPotential(int row,int column,double Potential)
{
    MatrixA.coeffRef(row,column)=-expl(Potential);

}
void Simul::AnalyzeMatrixPattern(void)
{
    MSolver.analyzePattern(MatrixA);
}
void Simul::FactorizeMatrix(void)
{
    MSolver.factorize(MatrixA);
}
void Simul::PrintConcen(char *fichero)
{
    FILE *pfile;
    pfile=fopen(fichero,"a");
    int kk,jj;
    for(kk=0;kk<MyData.NPointsGridSpecies;kk++)
    {
        for(jj=0;jj<MyData.NSpecies;jj++)
             fprintf(pfile,"%lf\t",Concentration[kk+jj*MyData.NPointsGridSpecies]);
        fprintf(pfile,"\n");
    }
    fprintf(pfile,"\n");
    fclose(pfile);
}

void Simul::SurfaceConditions(double Potential)
{
    int Ox,Red,Row,NumberOfSpecies;
    int SpeciesIndex[MAX_SPECIES];
    Ox=0;//species A
    Red=1;//species B
    switch(MyData.Mechanism)
    {
    case EMECHANISM:
        Ox=0;//species A
        Red=1;//species B
        //Nernst condition in row '0'
        Row=0;
        Nernst(Row,Ox,Red,Potential);

        //equal flux of Ox and Red in row 'NPointsGridSpecies'
        Row=MyData.NPointsGridSpecies;
        NumberOfSpecies=2;
        SpeciesIndex[0]=Ox;
        SpeciesIndex[1]=Red;
        Fluxes(Row,NumberOfSpecies,SpeciesIndex);
        break;
    case ECMECHANISM:
        Ox=0;//species A
        Red=1;//species B
        //Nernst condition in row '0'
        Row=0;
        Nernst(Row,Ox,Red,Potential);

        //equal flux of Ox and Red in row 'NPointsGridSpecies'
        Row=MyData.NPointsGridSpecies;
        NumberOfSpecies=2;
        SpeciesIndex[0]=Ox;
        SpeciesIndex[1]=Red;
        Fluxes(Row,NumberOfSpecies,SpeciesIndex);
        //zero flux of C (species 2) in row '2*NPointsGridSpecies'
        Row=2*MyData.NPointsGridSpecies;
        NumberOfSpecies=1;
        SpeciesIndex[0]=2; //species C
        Fluxes(Row,NumberOfSpecies,SpeciesIndex);
        break;
        //other cases can be added

    }


}
void Simul::Nernst(int Row,int Ox,int Red,double Potential)
{
    MContainer.push_back(MTriplet(Row,Ox*MyData.NPointsGridSpecies,1.));
    MContainer.push_back(MTriplet(Row,Red*MyData.NPointsGridSpecies,-expl(Potential)));

}
void Simul::Fluxes(int Row,int NEspecies,int *Index)
{
    int jj,kk;
    for(jj=0;jj<NEspecies;jj++)
    {
        for (kk=0;kk<MyData.PointsCurrent;kk++)
                MContainer.push_back(MTriplet(Row,Index[jj]*MyData.NPointsGridSpecies+kk,CoefBeta0[kk]/MyData.hh));
    }
}
