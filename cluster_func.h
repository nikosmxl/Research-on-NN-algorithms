#include <iostream>
#include <vector>
#include <bits/stdc++.h>

template<typename T>
std::string print_vector(std::vector<T> vect){
    std::stringstream temp;
    for (int i = 0 ; i < vect.size() ; i++){
        temp << vect[i];
        if (i != vect.size() - 1){
            temp << ", ";
        }
    }
    return temp.str();
}

template<typename T>
class Cluster{
    private:
        T** pixels;
        T* center;
        int DIM;
        std::size_t size;
    public:
        Cluster(T** p, int DIMENSION, int cntr) : pixels(p), size(0), DIM(DIMENSION) {
            center = new int[DIM];
            
            for (int i = 0 ; i < DIM ; i++){
                center[i] = pixels[cntr][i];
            }
        }

        ~Cluster(){
            delete[] center;
        }

        void update_center(int* sums){
            for (int i = 0 ; i < DIM ; i++){
                center[i] = sums[i] / size;
            }
        }

        T* get_center() const {
            return center;
        }

        void size_up(){
            size++;
        }

        void size_down(){
            size--;
        }

        int get_size() const{
            return size;
        }

        std::string print_center_coordinates() const{
            std::stringstream coordinates;
            for (int i = 0 ; i < DIM ; i++){
                coordinates << center[i];
                if (i != DIM - 1){
                    coordinates << ", ";
                }
            }
            return coordinates.str();
        }
        
};

template<typename T>
float R_init(Cluster<T>** clusters, int number_of_clusters, int DIMENSION, double (*dist)(T*, T*, int, int)){
    float min = dist(clusters[0]->get_center(),clusters[1]->get_center(),2,DIMENSION);
    for(int i = 0 ; i < number_of_clusters ; i++){
        for (int j = 0 ; j < number_of_clusters ; j++){
            if (i == j || (i == 0 && j == 1)) continue;

            double d = dist(clusters[i]->get_center(),clusters[j]->get_center(),2,DIMENSION);
            if (d < min){
                min = d;
            }
        }
    }
    return min;
}

template<typename T>
void assignment_lloyds(T** pixels, Cluster<T>** clusters, int* clustindexforps, int number_of_clusters, int NO_IMAGES, int DIMENSION, double (*dist)(T*, T*, int, int)){
    for (int i = 0 ; i < NO_IMAGES ; i++){
        int min = dist(pixels[i], clusters[0]->get_center(), 2, DIMENSION);
        int mincenter = 0;
        for(int j = 1; j < number_of_clusters; j++){
            double ddd = dist(pixels[i], clusters[j]->get_center(), 2, DIMENSION);
            if( ddd < min ){
                min = ddd;
                mincenter = j;
            }
        }
        clustindexforps[i] = mincenter;
        clusters[clustindexforps[i]]->size_up();    // Αυξανουμε το size του cluster που πηρε το vector
    }
}

template<typename T>
void kmeans_plusplus(T** pixels, Cluster<T>** clusters, int* centers, int number_of_clusters, int NO_IMAGES, int DIMENSION, double (*dist)(T*, T*, int, int)){
    centers[0] = rand() % NO_IMAGES;
    clusters[0] = new Cluster<T>(pixels, DIMENSION, centers[0]);

    double distances[NO_IMAGES];
    for (int i = 0 ; i < NO_IMAGES ; i++){
        distances[i] = dist(pixels[i], clusters[0]->get_center(), 2, DIMENSION);
    }

    for (int i = 1 ; i < number_of_clusters ; i++){
        double prob[NO_IMAGES];
        long double sum_square = 0;
        for (int j = 0 ; j < NO_IMAGES ; j++){
            bool already_center = false;
            for (int k = 0 ; k < number_of_clusters ; k++){
                if (centers[k] == -1) break;
                if (i == centers[k]){
                    already_center = true;
                    break;
                }
            }
            if (already_center) continue;

            prob[j] = pow(distances[j], 2);    // Di^2
            sum_square += prob[j];
        }

        for (int j = 0 ; j < NO_IMAGES ; j++){
            bool already_center = false;
            for (int k = 0 ; k < number_of_clusters ; k++){
                if (centers[k] == -1) break;
                if (i == centers[k]){
                    already_center = true;
                    break;
                }
            }
            if (already_center) continue;

            prob[j] /= sum_square;          // Επιλεγουμε το κεντρο του iοστου cluster αποδιδοντας πιθανοτητες σε καθε ενα συμφωνα με την kmeans++
        }

        std::default_random_engine generator;
        std::discrete_distribution<int> discrete(prob, prob + NO_IMAGES);

        int next_center = discrete(generator);  // Μεσω του discrete distribution θα επιλεξουμε ενα σημειο για κεντρο του επομενου cluster
        centers[i] = next_center;
        clusters[i] = new Cluster<T>(pixels, DIMENSION, centers[i]);

        if (i == number_of_clusters - 1) break;
        
        for (int j = 0 ; j < NO_IMAGES ; j++){
            double d = dist(pixels[j], clusters[i]->get_center(), 2, DIMENSION);
            if (d < distances[j]){
                distances[j] = d;         // Ενημερωνουμε τα distances αφου προστεθηκε ενα νεο cluster
            }
        }
    }
}

template<typename T>
void update(int number_of_clusters, int DIM, int picsnum, int* clustindexforps, T** vects, Cluster<T>** clusters){
    for(int i = 0; i < number_of_clusters; i++){
        int* sum = new int[DIM];
        for(int j = 0; j < DIM; j++){
            sum[j] = 0;
        }
        for(int j = 0; j < picsnum; j++){
            if( clustindexforps[j] == i ){
                for(int m = 0; m < DIM; m++){
                    sum[m] += vects[j][m];
                }
            }
        }
        clusters[i]->update_center(sum);
        delete[] sum;
    }
}

template<typename T>
void cluster_lloyds(T** pixels, Cluster<T>** clusters, int* clustindexforps, int number_of_clusters, int NO_IMAGES, int DIMENSION, float stopfactor, double (*dist)(T*, T*, int, int)){
    // Lloyd's algorithm
    // Assignment (Expectation)
    assignment_lloyds(pixels, clusters, clustindexforps, number_of_clusters, NO_IMAGES, DIMENSION, dist);

    // Update (Maximization)
    double changefactor = 100;
    while( changefactor >= stopfactor ){       // Οσο υπαρχει διαφορα πανω απο 5% σε σχεση με την τελευταια κατασταση συνεχιζουμε να ενημερωνουμε
        update(number_of_clusters, DIMENSION, NO_IMAGES, clustindexforps, pixels, clusters);
        int sumofchanged = 0;           // Αθροισμα διανυσματων που αλλαξαν cluster
        for(int i = 0; i < NO_IMAGES; i++){
            int min = dist(pixels[i], clusters[0]->get_center(), 2, DIMENSION);
            int mincenter = 0;
            for(int j = 1; j < number_of_clusters; j++){
                double ddd = dist(pixels[i], clusters[j]->get_center(), 2, DIMENSION);
                if( ddd < min ){
                    min = ddd;
                    mincenter = j;
                }
            }
            if( clustindexforps[i] != mincenter ){
                clusters[clustindexforps[i]]->size_down();  // Μειωνουμε το size του αφου θα χασει ενα vector
                clustindexforps[i] = mincenter;
                clusters[clustindexforps[i]]->size_up();    // Αυξανουμε το size του cluster που πηρε το vector
                sumofchanged++;
            }
        }
        changefactor = sumofchanged*100/(double)(NO_IMAGES);
    }
}

template<typename T>
void cluster_lsh(T** pixels, Cluster<T>** clusters, std::unordered_multimap<int, int>** mm, int* clustindexforps, int* centers, bool* assigned, int** w, double** t, int** rs, long* id, float stopfactor, int number_of_clusters, int L, int K, long M, int NO_IMAGES, int DIMENSION, float R, double (*dist)(T*, T*, int, int)){
    double changefactor = 100;
    while (changefactor >= stopfactor){        // Οσο υπαρχει διαφορα πανω απο stopfactor% σε σχεση με την τελευταια κατασταση συνεχιζουμε να ενημερωνουμε
        bool checked[NO_IMAGES];
        bool changed[NO_IMAGES];
        for (int i = 0 ; i < NO_IMAGES ; i++){
            changed[i] = false;
            checked[i] = false;
        }
        int sumofchanged = 0;           // Αθροισμα διανυσματων που αλλαξαν cluster
        int sumofchecked = 0;
        for(int l = 0; l < L; l++){
            for (int j = 0 ; j < number_of_clusters ; j++){
                int key = g(pixels, w[l], t[l], rs[l], id, centers[j], K, M, NO_IMAGES/8, DIMENSION, l);
                auto itr = mm[l]->equal_range(key);
                for (auto it = itr.first; it != itr.second; it++) {
                    double d = dist(pixels[it->second],clusters[j]->get_center(),2,DIMENSION);
                    if ( d >= R ) continue;

                    checked[it->second] = true;
                    if (assigned[it->second]){
                        if (d < dist(pixels[it->second], clusters[clustindexforps[it->second]]->get_center(), 2, DIMENSION)){
                            clusters[clustindexforps[it->second]]->size_down();
                            clustindexforps[it->second] = j;
                            clusters[clustindexforps[it->second]]->size_up();
                            changed[it->second] = true;
                        }
                    }
                    else{
                        clustindexforps[it->second] = j;
                        clusters[clustindexforps[it->second]]->size_up();
                        assigned[it->second] = true;
                        changed[it->second] = true;
                    }
                }
            }
        }
        R *= 2;
        for (int i = 0 ; i < NO_IMAGES ; i++){
            if (checked[i]) sumofchecked++;
            if (changed[i]) sumofchanged++;
        }
        changefactor = sumofchanged*100/(double)(sumofchecked);
    }
}

template<typename T>
void cluster_hypercube(T** pixels, Cluster<T>** clusters, std::unordered_multimap<long, int>& hypercube, long* query_key, int* clustindexforps, bool* assigned, int number_of_clusters, float stopfactor, int dt, int probes, int Mcube, float R, int NO_IMAGES, int DIMENSION, double (*dist)(T*, T*, int, int)){
    double changefactor = 100;
    while (changefactor >= stopfactor){        // Οσο υπαρχει διαφορα πανω απο stopfactor% σε σχεση με την τελευταια κατασταση συνεχιζουμε να ενημερωνουμε
        bool checked[NO_IMAGES];
        bool changed[NO_IMAGES];
        for (int i = 0 ; i < NO_IMAGES ; i++){
            changed[i] = false;
            checked[i] = false;
        }
        int sumofchanged = 0;           // Αθροισμα διανυσματων που αλλαξαν cluster
        int sumofchecked = 0;
        for (int i = 0 ; i < number_of_clusters ; i++){
            auto itr = hypercube.equal_range(query_key[i]);

            std::bitset<32> binary_query_key(query_key[i]); // Convert p to binary representation
            int countVerticesRange = 0;
            int countCubeElementsRange = 0;
            int hamming_dist_wantedRange = 0;
            int hd = 0;

            while (countCubeElementsRange < Mcube && countVerticesRange < probes){
                for (auto it = itr.first; it != itr.second; it++) {
                    countCubeElementsRange++;
                    double d = dist(pixels[it->second],clusters[i]->get_center(),2,DIMENSION);
                    if(d >= R) continue;
                    
                    checked[it->second] = true;
                    if (assigned[it->second]){
                        if (d < dist(pixels[it->second], clusters[clustindexforps[it->second]]->get_center(), 2, DIMENSION)){
                            clusters[clustindexforps[it->second]]->size_down();
                            clustindexforps[it->second] = i;
                            clusters[clustindexforps[it->second]]->size_up();
                            changed[it->second] = true;
                        }
                    }
                    else{
                        clustindexforps[it->second] = i;
                        clusters[clustindexforps[it->second]]->size_up();
                        assigned[it->second] = true;
                        changed[it->second] = true;
                    }
                }
                countVerticesRange++;

                if (hd == query_key[i]){
                    hd++;
                    hd = hd % dt;
                }
                if (hd == 0) hamming_dist_wantedRange++;
                while (1){
                    std::bitset<32> binary_hd(hd); // Convert num to binary representation
                    int differingBits = (binary_query_key ^ binary_hd).count(); // Count differing bits
                    if (differingBits == hamming_dist_wantedRange){
                        itr = hypercube.equal_range(hd);
                        break;
                    }

                    hd++;
                    hd = hd % dt;
                    if (hd == 0) hamming_dist_wantedRange++;
                }
                
            }
        }
        
        R *= 2;
        for (int i = 0 ; i < NO_IMAGES ; i++){
            if (checked[i]) sumofchecked++;
            if (changed[i]) sumofchanged++;
        }
        changefactor = sumofchanged*100/(double)(sumofchecked);
    }
}

template<typename T>
void assign_unassigned(T** pixels, Cluster<T>** clusters, int* clustindexforps, bool* assigned, int number_of_clusters, int NO_IMAGES, int DIMENSION, double (*dist)(T*, T*, int, int)){
    for (int i = 0 ; i < NO_IMAGES ; i++){
        if (assigned[i]) continue;

        int min = dist(pixels[i], clusters[0]->get_center(), 2, DIMENSION);
        int mincenter = 0;
        for(int j = 1; j < number_of_clusters; j++){
            double ddd = dist(pixels[i], clusters[j]->get_center(), 2, DIMENSION);
            if( ddd < min ){
                min = ddd;
                mincenter = j;
            }
        }

        clustindexforps[i] = mincenter;
        clusters[clustindexforps[i]]->size_up();
    }
}

int readconf(std::string , int& , int& , int& , int& , int& , int& );