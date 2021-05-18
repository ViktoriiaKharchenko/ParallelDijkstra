# include <time.h>
# include <omp.h>
#include<iostream>
using namespace std;

#define INFINITY 10000000
void Read_matrix(int mat[], int n);
void Dijkstra_Init(int mat[], int pred[], int dist[], int known[], int n);
void RandomDataInitialization(int* pMatrix, int n);
void Dijkstra(int mat[], int dist[], int pred[], int n);
int Find_min_dist(int start, int end, int dist[], int known[], int n, int * my_min);
void Print_matrix(int mat[], int rows, int cols);
void Print_dists(int dist[], int n);
void Print_paths(int pred[], int n);


int main()
{
    int* dist, * pred, * mat;
    int n;
    time_t  start, finish;
    double duration;
    cin >> n;
    mat = new int[n * n];
    dist = new int[n];
    pred = new int[n];
    //Read_matrix(mat, n);
    RandomDataInitialization(mat, n);
    start = clock();
    Dijkstra(mat, dist, pred, n);
    finish = clock();
    duration = (finish - start) / double(CLOCKS_PER_SEC);
    cout << "\n Time of execution: " << duration;
    //Print_dists(dist, n);
    //Print_paths(pred, n);

    delete[] mat;
    delete[] dist;
    delete[] pred;
}
void Read_matrix(int mat[], int n) {

    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            cin >> mat[i * n + j];
            if (mat[i * n + j] == 0) {
                mat[i * n + j] = INFINITY;
            }
        }
}

void RandomDataInitialization(int* pMatrix, int n) {
    int i, j, k;
    srand(unsigned(clock()));
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            k = rand() / (int)(100);
            if (k == 0) k = INFINITY;
            pMatrix[i * n + j] = k;
        }
    }
}
void Dijkstra_Init(int mat[], int pred[], int dist[], int known[], int n) {

#pragma omp parallel for 
    for (int i = 0; i < n; i++) {
        known[i] = 0;
        dist[i] = mat[i];
        pred[i] = 0;
    }
    
}
int Find_min_dist(int start, int end, int dist[], int known[], int n, int * my_min) {
    int shortest_dist = INFINITY;
    int u = -1;
    for (int i = start; i <= end; i++) {
        if (!known[i]) {
            if (dist[i] < shortest_dist) {
                shortest_dist = dist[i];
                u = i;
            }
        }
    }
    *my_min = shortest_dist;
    return u;
}
void Dijkstra(int mat[], int dist[], int pred[], int n) {

    int  v, u, new_dist;
    int* known = new int[n];
    int my_min;
    int glob_min = INFINITY;
    int glob_u = 0;
    int my_first;
    int proc_id;
    int my_last;
    int proc_num;


    Dijkstra_Init(mat, pred, dist, known, n);

#pragma omp parallel private ( my_first, proc_id, my_last, my_min, new_dist, u, v) \
  shared ( known, dist,pred, glob_min, glob_u, proc_num, mat )
    {
        proc_id = omp_get_thread_num();
        proc_num = omp_get_num_threads();
        my_first = (proc_id * n) / proc_num;
        my_last = ((proc_id + 1) * n) / proc_num - 1;

        for (v = 0; v < n - 1; v++)
        {
            u = Find_min_dist(my_first, my_last, dist, known, n, &my_min);

        #pragma omp critical
            {
                if (my_min < glob_min)
                {
                    glob_min = my_min;
                    glob_u = u;
                }
            }
        #pragma omp barrier

                for (int i = my_first; i <= my_last; i++) {
                    if (!known[i]) {
                        new_dist = glob_min + mat[glob_u * n + i];
                        if (new_dist < dist[i]) {
                            dist[i] = new_dist;
                            pred[i] = glob_u;
                        }
                    }
                }
        #pragma omp barrier

         }

   }
        delete[]known;
  }


void Print_matrix(int mat[], int rows, int cols) {
    int i, j;

    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++)
            if (mat[i * cols + j] == INFINITY)
                cout << "i ";
            else
                cout << mat[i * cols + j] << " ";
        cout << "\n";
    }

    cout << "\n";
}
void Print_dists(int dist[], int n) {
    int v;

    cout << "  v    dist 0->v\n";
    cout << "----   ---------\n";

    for (v = 1; v < n; v++) {
        if (dist[v] == INFINITY) {
            cout << v << "       inf\n";
        }
        else
            cout << v << "       " << dist[v] << "\n";
    }
    cout << "\n";
}


void Print_paths(int pred[], int n) {
    int v, w, * path, count, i;

    path = new int[n];

    cout << "  v     Path 0->v\n";
    cout << "----    ---------\n";
    for (v = 1; v < n; v++) {
        cout << v << ":    ";
        count = 0;
        w = v;
        while (w != 0) {
            path[count] = w;
            count++;
            w = pred[w];
        }
        cout << "0 ";
        for (i = count - 1; i >= 0; i--)
            cout << path[i] << " ";
        cout << "\n";
    }
    delete[]path;
}