#include "bp.h"

using namespace std;

/* Constructor */
bp::bp (int a, int b, int c, int d){
    k = a;
    m1 = b;
    n = c;
    m2 = d;


    if(n==0) // Single Bimodal Predictor
    {  
        bimodalRows = ceil(pow(2,m2));
        bimodalPredTable = new int[bimodalRows];
        bimodalMask = ((1u << (m2 + 2)) - 1) & ~((1u << 2) - 1);
    }
    else if (n!=0)
    {
        if(k==0) // gshare predictor
        {
            gshareRows = ceil(pow(2,m1));
            gsharePredTable = new int[gshareRows];
            gshareMask = (1u << (m1 - n)) - 1 ; 
            hRegMask = (1u << n) - 1;
        }
        else if(k!=0)
        {
            chooserRows = ceil(pow(2,k));
            chooserTable = new int[chooserRows];
            chooserMask = ((1u << (k + 2)) - 1) & ~((1u << 2) - 1);
            bimodalRows = ceil(pow(2,m2));
            bimodalPredTable = new int[bimodalRows];
            bimodalMask = ((1u << (m2 + 2)) - 1) & ~((1u << 2) - 1);
            gshareRows = ceil(pow(2,m1));
            gsharePredTable = new int[gshareRows];
            gshareMask = (1u << (m1 - n)) - 1 ; 
            hRegMask = (1u << n) - 1;
        }
    }
    
    predictions = 0;
    mispredictions = 0;
    initbp();
}

void bp::initbp()
{
    if(n==0)
    {
        for(int i=0;i<bimodalRows;i++)
        {
            bimodalPredTable[i] = 2;
        }
    }
    else if(n!=0)
    {
        if(k==0)
        {
            for(int i=0; i<gshareRows; i++)
            {
                gsharePredTable[i] = 2;
            }
            histReg = 0;
        }
        if(k!=0)
        {
            for(int i=0;i<bimodalRows;i++)
            {
                bimodalPredTable[i] = 2;
            }
            for(int i=0; i<gshareRows; i++)
            {
                gsharePredTable[i] = 2;
            }
            histReg = 0;
            for(int i=0; i<chooserRows;i++)
            {
                chooserTable[i] = 1;
            }
        }
    }
}

void bp::bimodalPrint()
{
    cout<<"FINAL BIMODAL CONTENTS"<< endl;
    for(int i=0; i<bimodalRows; i++)
    {
        cout<<" "<<setw(6)<<left<<i<<left<<bimodalPredTable[i]<<endl;
    }
}

void bp::gsharePrint()
{
    cout<<"FINAL GSHARE CONTENTS"<< endl;
    for(int i=0; i<gshareRows; i++)
    {
        cout<<" "<<setw(6)<<left<<i<<left<<gsharePredTable[i]<<endl;
    }

}

void bp::hybridPrint()
{
    cout<<"FINAL CHOOSER CONTENTS"<< endl;
    for(int i=0; i<chooserRows; i++)
    {
        cout<<" "<<setw(6)<<left<<i<<left<<chooserTable[i]<<endl;
    }
    cout<<"FINAL GSHARE CONTENTS"<< endl;
    for(int i=0; i<gshareRows; i++)
    {
        cout<<" "<<setw(6)<<left<<i<<left<<gsharePredTable[i]<<endl;
    }
    cout<<"FINAL BIMODAL CONTENTS"<< endl;
    for(int i=0; i<bimodalRows; i++)
    {
        cout<<" "<<setw(6)<<left<<i<<left<<bimodalPredTable[i]<<endl;
    }

}

void bp::predict(unsigned int PC, const char *tn)
{
    unsigned int a = 0;
    predictions++;
    int bimodalPrediction=-1; // bimodal Predicted Counter
    int gsharePrediction=-1;  // gshare Predicted Counter
    int overallPrediction=-1; // Prediction from the chooser 
    int bimodalResult=-1;
    int gshareResult=-1;
    int actual=-1; // Updated Counter
    PCaddr = PC;
    if(*tn == 't')
        {actual = 1;}    // actual taken
    else if(*tn == 'n')
        {actual = 0;}      // actual not taken
  
    if(n==0) // Single Bimodal Counter
    {
        index = (PCaddr & bimodalMask) >> 2;
        if(bimodalPredTable[index] >= 2)
            { bimodalPrediction = 1;}                            // Predict taken
        else if (bimodalPredTable[index] < 2)
            { bimodalPrediction = 0;}                            // Predict not-taken

        if (bimodalPrediction != actual)
        { mispredictions++;}

        // Updating the counter 
        if (actual == 1)
        {
            if(bimodalPredTable[index] < 3)
            {
                bimodalPredTable[index] = bimodalPredTable[index] + 1;
            }
        }
        else if (actual == 0)
        {
            if(bimodalPredTable[index] != 0)
            {
                bimodalPredTable[index] = bimodalPredTable[index] - 1;
            }
        }
    }
    else if (n!=0)
    {
        if(k==0) // gshare counter  
        {
            a = histReg ^ ( (PCaddr >> ( 2 + m1 - n)) & hRegMask );
	        index = (a << (m1 - n)) + ( (PCaddr >> 2) & gshareMask );
            if(gsharePredTable[index] >= 2)
                { gsharePrediction = 1;}                            // Predict taken
            else if (gsharePredTable[index] < 2)
                { gsharePrediction = 0;}                            // Predict not-taken
            
            if (gsharePrediction != actual)
                { mispredictions++;}
            
            // Updating the counter 
            if (actual == 1)
            {
                if(gsharePredTable[index] < 3)
                {
                    gsharePredTable[index] = gsharePredTable[index] + 1;
                }
            }
            else if (actual == 0)
            {
                if(gsharePredTable[index] != 0)
                {
                    gsharePredTable[index] = gsharePredTable[index] - 1;
                }
            }
            histReg = (actual << (n-1))  | (histReg >> 1);
        }

        else if (k!=0) // Hybrid Predictor
        {
            // Bimodal Prediction 
            index = (PCaddr & bimodalMask) >> 2;
            if(bimodalPredTable[index] >= 2)
                { bimodalPrediction = 1;}                            // Predict taken
            else if (bimodalPredTable[index] < 2)
                { bimodalPrediction = 0;}                            // Predict not taken

            if ( bimodalPrediction == actual)                        // Getting bimodal predictor result
		        {bimodalResult = 1;}
	        else {bimodalResult = 0;}    

            //ghsare Prediction
            a = histReg ^ ( (PCaddr >> ( 2 + m1 - n)) & hRegMask );
	        index = (a << (m1 - n)) + ( (PCaddr >> 2) & gshareMask );
            if(gsharePredTable[index] >= 2)
                { gsharePrediction = 1;}                            // Predict taken
            else if (gsharePredTable[index] < 2)
                { gsharePrediction = 0;}                            // Predict not-taken

            if ( gsharePrediction == actual)                        // Getting gshare predictor result
		        {gshareResult = 1;}
	        else {gshareResult = 0;}  

            //Overall Prediction
            index = (PCaddr & chooserMask) >> 2;
            if(chooserTable[index] >= 2)
                { overallPrediction = gsharePrediction;}                            // Choose gshare
            else if (chooserTable[index] < 2)
                { overallPrediction = bimodalPrediction;}                            // Choose bimodal

            if (overallPrediction != actual)
                { mispredictions++;}    


            //Updating 
            if(chooserTable[index] >= 2) // gshare
                { 
                    index = (a << (m1 - n)) + ( (PCaddr >> 2) & gshareMask );
                    if (actual == 1)
                    {
                        if(gsharePredTable[index] < 3)
                        {
                            gsharePredTable[index] = gsharePredTable[index] + 1;
                        } 
                    }
                    else if (actual == 0)
                    {
                        if(gsharePredTable[index] != 0)
                        {
                            gsharePredTable[index] = gsharePredTable[index] - 1;
                        }
                    }

                }                            
            else if (chooserTable[index] < 2) // bimodal
                {
                    index = (PCaddr & bimodalMask) >> 2;
                    if (actual == 1)
                    {
                        if(bimodalPredTable[index] < 3)
                        {
                            bimodalPredTable[index] = bimodalPredTable[index] + 1;
                        }
                    }
                    else if (actual == 0)
                    {
                        if(bimodalPredTable[index] != 0)
                        {
                            bimodalPredTable[index] = bimodalPredTable[index] - 1;
                        }
                    }
                }  
            histReg = (actual << (n-1))  | (histReg >> 1); // Global Branch History Register

            // Chooser table 
            index = (PCaddr & chooserMask) >> 2;
            if (( gshareResult == 1 ) & ( bimodalResult == 0))
            {
                if(chooserTable[index] < 3)
                {
                    chooserTable[index] = chooserTable[index] + 1;
                }
            }
            else if (( gshareResult == 0 ) & ( bimodalResult == 1))
            {
                if(chooserTable[index] != 0)
                {
                    chooserTable[index] = chooserTable[index] - 1;
                }
            }

        }
        
    }
     
}