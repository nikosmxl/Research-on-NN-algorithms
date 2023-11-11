#include "cluster_func.h"

std::string print_vector(std::vector<int> vect){
    std::stringstream temp;
    for (int i = 0 ; i < vect.size() ; i++){
        temp << vect[i];
        if (i != vect.size() - 1){
            temp << ", ";
        }
    }
    return temp.str();
}

Cluster::Cluster(int** p, int DIMENSION, int cntr) : pixels(p), size(0), DIM(DIMENSION) {
    center = new int[DIM];
    
    for (int i = 0 ; i < DIM ; i++){
        center[i] = pixels[cntr][i];
    }
}

Cluster::~Cluster(){
    delete[] center;
}

void Cluster::update_center(int* sums){
    for (int i = 0 ; i < DIM ; i++){
        center[i] = sums[i] / size;
    }
}

int* Cluster::get_center() const {
    return center;
}

void Cluster::size_up(){
    size++;
}

void Cluster::size_down(){
    size--;
}

int Cluster::get_size() const{
    return size;
}

std::string Cluster::print_center_coordinates() const{
    std::stringstream coordinates;
    for (int i = 0 ; i < DIM ; i++){
        coordinates << center[i];
        if (i != DIM - 1){
            coordinates << ", ";
        }
    }
    return coordinates.str();
}
void assignment_lloyds(int** pixels, Cluster** clusters, int* clustindexforps, int number_of_clusters, int NO_IMAGES, int DIMENSION, double (*dist)(int*, int*, int, int)){
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

float R_init(Cluster** clusters, int number_of_clusters, int DIMENSION, double (*dist)(int*, int*, int, int)){
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

void kmeans_plusplus(int** pixels, Cluster** clusters, int* centers, int number_of_clusters, int NO_IMAGES, int DIMENSION, double (*dist)(int*, int*, int, int)){
    centers[0] = rand() % NO_IMAGES;
    clusters[0] = new Cluster(pixels, DIMENSION, centers[0]);

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
        clusters[i] = new Cluster(pixels, DIMENSION, centers[i]);

        if (i == number_of_clusters - 1) break;
        
        for (int j = 0 ; j < NO_IMAGES ; j++){
            double d = dist(pixels[j], clusters[i]->get_center(), 2, DIMENSION);
            if (d < distances[j]){
                distances[j] = d;         // Ενημερωνουμε τα distances αφου προστεθηκε ενα νεο cluster
            }
        }
    }
}

void update(int number_of_clusters, int DIM, int picsnum, int* clustindexforps, int** vects, Cluster** clusters){
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

int readconf(std::string conf, int& number_of_clusters, int& L, int& K, int& Mcube, int& dt, int& probes){
    std::ifstream confFile(conf); // Open the file for reading
    if (!confFile) {
        std::cerr << "Error: Could not open the file." << std::endl;
        return 1;
    }

    std::string line;
    int i = 0;
    while (std::getline(confFile, line)) { // Read each line from the file
        char* token = std::strtok(&line[0], " \t\n"); // Tokenize the line using space, tab, and newline as delimiters
        int flag = -1;
        while (token != nullptr) {
            if (strcmp(token, "//") == 0) break;

            if (flag == -1){
                if (strcmp(token, "number_of_clusters:") == 0){
                    flag = 0;
                }
                if (strcmp(token, "number_of_vector_hash_tables:") == 0){
                    flag = 1;
                }
                if (strcmp(token, "number_of_vector_hash_functions:") == 0){
                    flag = 2;
                }
                if (strcmp(token, "max_number_M_hypercube:") == 0){
                    flag = 3;
                }
                if (strcmp(token, "number_of_hypercube_dimensions:") == 0){
                    flag = 4;
                }
                if (strcmp(token, "number_of_probes:") == 0){
                    flag = 5;
                }
            }
            else{
                if (flag == 0){
                    number_of_clusters = atoi(token);
                }
                else if (flag == 1){
                    L = atoi(token);
                }
                else if (flag == 2){
                    K = atoi(token);
                }
                else if (flag == 3){
                    Mcube = atoi(token);
                }
                else if (flag == 4){
                    dt = atoi(token);
                }
                else if (flag == 5){
                    probes = atoi(token);
                }
                flag = -1;
            }
            token = std::strtok(nullptr, " \t\n"); // Get the next token
        }
        i++;
    }

    confFile.close(); // Close the file
    return 0;
}