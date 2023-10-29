#include <iostream>
#include <fstream>
#include <bits/stdc++.h>

int readfile(std::string, int** , int , int );

double dist(int* ,int* , int , int );

class Neibs
{
    private:
        int** p;
        int** queries;
        int DIM;
        int* nn;
        int size;
        int maxsize;
        int q;
        double (*distfunc)(int*, int*, int, int);
    public:
    Neibs(int** pixels, int** qs, int dimension, int s,int query, double (*func)(int*, int*, int, int));

    ~Neibs();

    int* givenn() const;

    int givenn(int i) const;

    int give_size() const;

    void insertionsort_insert(int v);

    double givedist(int i) const;
};