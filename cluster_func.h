#include <iostream>
#include <vector>
#include <bits/stdc++.h>

std::string print_vector(std::vector<int> );

class Cluster{
    private:
        int** pixels;
        int* center;
        int DIM;
        int size;
    public:
        Cluster(int** p, int DIMENSION, int cntr);

        ~Cluster();

        void update_center(int* sums);

        int* get_center() const ;

        void size_up();

        void size_down();

        int get_size() const;

        std::string print_center_coordinates() const;
        
};

float R_init(Cluster** , int , int , double (*dist)(int*, int*, int, int));

void assignment_lloyds(int** , Cluster** , int* , int , int , int , double (*dist)(int*, int*, int, int));

void kmeans_plusplus(int** , Cluster** , int* , int , int , int , double (*dist)(int*, int*, int, int));

void update(int , int , int , int* , int** , Cluster** );

int readconf(std::string , int& , int& , int& , int& , int& , int& );