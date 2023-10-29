#include <iostream>
#include <random>
#include "generals.h"

int main(){
    //floyd kmeans
    const int NO_QUERIES = 10;
    const int NO_IMAGES = 60000;
    const int DIMENSION = 784;

    int numofcenters = 4;
    int** centers = new int*[numofcenters];
    int** maxNminD = new int*[NO_IMAGES];
    int* clustindexforps = new int[NO_IMAGES];

    int** pixels = new int*[NO_IMAGES];
    int** queries = new int*[NO_QUERIES];

    for (int i = 0 ; i < NO_IMAGES ; i++){
        pixels[i] = new int[DIMENSION];
    }

    for (int i = 0 ; i < NO_QUERIES ; i++){
        queries[i] = new int[DIMENSION];
    }

    // maxNminD exi gia kathe diastasi tin max ke ti min ti apo ola ta dianimsata
    for(int i = 0; i < DIMENSION; i++){ //0 max    1 min
        maxNminD[i] = new int[2];
        maxNminD[i][0] = pixels[0][i];
        maxNminD[i][1] = pixels[0][i];
        for(int j = 1; j <NO_IMAGES; j++){
            if( pixels[j][i] > maxNminD[i][0] ){
                maxNminD[i][0] = pixels[j][i];
            }
            if( pixels[j][i] < maxNminD[i][1] ){
                maxNminD[i][1] = pixels[j][i];
            }
        }
    }
    // to kathe cluster tha exi tixees times ti kathe diastasi apo min eos max aftis tis diastasis olon ton dianismaton
    for(int i = 0; i < numofcenters; i++){
        centers[i] = new int[DIMENSION];
        for(int j = 0; j < DIMENSION; j++){
            int modder = maxNminD[i][0]-maxNminD[i][1]+1;   //max-min+1
            centers[i][j] = rand()%modder + maxNminD[i][1]; //rand()%modder + min  ara tixea timi apo min eos max
        }
    }


    /*for(int i = 0; i < picsnum+1; i++){
        int min = dist(vects[i],centers[0]);
        int mincenter = 0;
        for(int j = 1; j < numofcenters; j++){
            if( int ddd = dist(vects[i],centers[j]) < min ){
                min = ddd;
                mincenter = j;
            }
        }
        clustindexforps[i] = mincenter;
        //gia to i dianisma vazume os cluster tu to mincenter
    }*/

    //ksanavriskoyme ta kedra me meso oro

    //update();
    // proto mirasma sta centers
    for(int i = 0; i < NO_IMAGES; i++){
        double min = dist(pixels[i],centers[0]);
        int mincenter = 0;
        for(int j = 1; j < numofcenters; j++){
            double ddd = dist(pixels[i],centers[j]);
            if( ddd < min ){
                min = ddd;
                mincenter = j;
            }
         }
            
        clustindexforps[i] = mincenter;
                
    }


    // oso yparxei pano apo 5% diafora apo prin sinexizoyme na kanume update
    double changefactor = 100;

    while( changefactor >= 5 ){
        update();
        int sumofchanged = 0;
        for(int i = 0; i < NO_IMAGES; i++){
            int min = dist(pixels[i],centers[0]);
            int mincenter = 0;
            for(int j = 1; j < numofcenters; j++){
                double ddd = dist(pixels[i],centers[j]);
                if( ddd < min ){
                    min = ddd;
                    mincenter = j;
                }
            }
            if( clustindexforps[i] != mincenter ){
                clustindexforps[i] = mincenter;
                sumofchanged++;
            }
        }
        changefactor = sumofchanged*100/(NO_IMAGES);
    }

    //floyd kmeans
}



void update(int numofcenters,int DIM,int picsnum,int* clustindexforps,int** vects,int* centers){
    for(int i = 0; i < numofcenters; i++){
        int* sum = new int[DIM];
        for(int j = 0; j < DIM; j++){
            sum[j] = 0;
        }
        int clustersize = 0;
        for(int j = 0; j < picsnum+1; j++){
            if( clustindexforps[j] == i ){
                clustersize++;
                for(int m = 0; m < DIM; m++){
                    sum[m] += vects[j][m];
                }
            }
        }
        if( clustersize > 0 ){
            for(int j = 0; j < vectdim(q); j++){
                sum[j] = sum[j]/clustersize++;
            }
            centers[i] = sum;
        }
        
    }
}