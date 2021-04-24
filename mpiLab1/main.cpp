#include <time.h>
#include<iostream>
#include <mpi.h>
#define INFINITY 10000000
using namespace std;

int Read_n(int my_rank, MPI_Comm comm, int proc_num, int *loc_n, int* send_ind);
void Read_matrix(int loc_mat[], int n, int loc_n,
    int my_rank, MPI_Comm comm, int proc_num, int * send_ind, int * dim, double * start);
void RandomDataInitialization(int* pMatrix, int n);

void Dijkstra_Init(int loc_mat[], int loc_pred[], int loc_dist[], int loc_known[],
    int my_rank, int loc_n);
void Dijkstra(int loc_mat[], int loc_dist[], int loc_pred[], int loc_n, int n,
    MPI_Comm comm, int send_ind);
int Find_min_dist(int loc_dist[], int loc_known[], int loc_n);
void Print_matrix(int global_mat[], int rows, int cols);
void Print_dists(int global_dist[], int n);
void Print_paths(int global_pred[], int n);

int main(int argc, char** argv) {
    int* loc_mat, * loc_dist, * loc_pred, * global_dist = NULL, * global_pred = NULL, *recv_num, *displ;
    int my_rank, p, loc_n, n, send_ind, loc_dim, loc_dim_ind;
    MPI_Comm comm;
    double Start, Finish, Duration;
    MPI_Init(NULL, NULL);
    comm = MPI_COMM_WORLD;
    MPI_Comm_rank(comm, &my_rank);
    MPI_Comm_size(comm, &p);
    recv_num = new int[p];
    displ = new int[p];

    n = Read_n(my_rank, comm, p, &loc_n, &send_ind);
    loc_dim = n * loc_n;
    loc_dim_ind = send_ind * n;
   // cout << "Rank = " << my_rank << " col_num = " << loc_n << endl;

    loc_mat = new int[n * loc_n];
    loc_dist = new int[loc_n];
    loc_pred = new int[loc_n];

    if (my_rank == 0) {
        global_dist = new int[n];
        global_pred = new int[n];
    }
    Read_matrix(loc_mat, n, loc_n, my_rank, comm,p, &loc_dim_ind, &loc_dim, &Start);

    //Start = MPI_Wtime();

    Dijkstra(loc_mat, loc_dist, loc_pred, loc_n, n, comm, send_ind);

    MPI_Allgather(&loc_n, 1, MPI_INT, recv_num, 1, MPI_INT, comm);
    MPI_Allgather(&send_ind, 1, MPI_INT, displ, 1, MPI_INT, comm);

    MPI_Gatherv(loc_dist, loc_n, MPI_INT, global_dist,recv_num, displ, MPI_INT, 0, comm);
    MPI_Gatherv(loc_pred, loc_n, MPI_INT, global_pred, recv_num, displ, MPI_INT, 0, comm);
         
    Finish = MPI_Wtime();
    Duration = Finish - Start;

    if (my_rank == 0) {
        //Print_dists(global_dist, n);
       // Print_paths(global_pred, n);
        cout << "\nTime of execution = "<< Duration << endl;
        delete[] global_dist;
        delete[] global_pred;
    }

    delete[] loc_mat;
    delete[] loc_pred;
    delete[] loc_dist;
    MPI_Finalize();
    return 0;
}

void RandomDataInitialization(int* pMatrix, int n) {
    int i, j,k; 
    srand(unsigned(clock()));
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            k = rand() / (int)(100);
            if (k == 0) k = INFINITY;
            pMatrix[i * n + j] = k;
        }
           
    }
}
int Read_n(int my_rank, MPI_Comm comm, int proc_num, int* loc_n, int* send_ind) {
    int n, *col_num = NULL, *ind_num = NULL;

    if (my_rank == 0)
        cin >> n;

    MPI_Bcast(&n, 1, MPI_INT, 0, comm);
    

    if (my_rank == 0) {
        col_num = new int[proc_num];
        ind_num = new int[proc_num];

        int RestRows = n % proc_num;
        for (int i = 0; i < proc_num; i++) {
            col_num[i] = n / proc_num ;
            if (RestRows != 0) {
                col_num[i]++;
                RestRows--;
            }
        }
        ind_num[0] = 0;
        for (int i = 1; i < proc_num; i++)
            ind_num[i] = ind_num[i - 1] + col_num[i - 1];
    }
    MPI_Scatter(col_num, 1, MPI_INT, loc_n, 1, MPI_INT, 0, comm);
    MPI_Scatter(ind_num, 1, MPI_INT, send_ind, 1, MPI_INT, 0, comm);

    if (my_rank == 0) {
        delete[]col_num;
        delete[]ind_num;
    }
    return n;
}

void Read_matrix(int loc_mat[], int n, int loc_n, int my_rank, MPI_Comm comm, int proc_num, int* send_int, int* dim, double * Start) {
    int* mat = NULL, i, j, * mat_col = NULL;
    int* col_num = new int[proc_num];
    int* ind_num = new int[proc_num];
    int* loc_mat2 = new int[n * loc_n];

   /* if (my_rank == 0) {
        mat = new int[n * n];
        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++) {
                cin >> mat[i * n + j];
                if (mat[i * n + j] == 0) {
                    mat[i * n + j] = INFINITY;
                }
            }
        mat_col = new int[n * n];
        for (i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                mat_col[j * n + i] = mat[i * n + j];
            }
        }
    }*/


    if (my_rank == 0) {
        mat_col = new int[n * n];
        RandomDataInitialization(mat_col, n);
    }

    *Start = MPI_Wtime();

    MPI_Allgather(dim, 1, MPI_INT, col_num, 1, MPI_INT, comm);
    MPI_Allgather(send_int, 1, MPI_INT, ind_num, 1, MPI_INT, comm);


    MPI_Scatterv(mat_col, col_num, ind_num, MPI_INT, loc_mat2, n * loc_n, MPI_INT, 0, comm);

            for (int i = 0; i < n; i++) {
                for (int j = 0; j < loc_n; j++) {
                    loc_mat[i * loc_n + j] = loc_mat2[j * n + i];
                }
            }
        
    

    if (my_rank == 0) {
        delete[]mat;
        delete[]mat_col;

    }
    delete[]col_num;
    delete[]ind_num;
    delete[]loc_mat2;
}

void Dijkstra_Init(int loc_mat[], int loc_pred[], int loc_dist[], int loc_known[],
    int my_rank, int loc_n) {
    int loc_v;

    if (my_rank == 0)
        loc_known[0] = 1;
    else
        loc_known[0] = 0;

    for (loc_v = 1; loc_v < loc_n; loc_v++)
        loc_known[loc_v] = 0;

    for (loc_v = 0; loc_v < loc_n; loc_v++) {
        loc_dist[loc_v] = loc_mat[0 * loc_n + loc_v];
        loc_pred[loc_v] = 0;
    }
}

void Dijkstra(int loc_mat[], int loc_dist[], int loc_pred[], int loc_n, int n,
    MPI_Comm comm, int ind_num) {

    int i, loc_v, loc_u, glbl_u, new_dist, my_rank, dist_glbl_u;
    int* loc_known;
    int my_min[2];
    int glbl_min[2];

    MPI_Comm_rank(comm, &my_rank);
    loc_known = new int[loc_n];

    Dijkstra_Init(loc_mat, loc_pred, loc_dist, loc_known, my_rank, loc_n);

   
    for (i = 0; i < n - 1; i++) {
       
        loc_u = Find_min_dist(loc_dist, loc_known, loc_n);
        if (loc_u != -1) {
            my_min[0] = loc_dist[loc_u];
            my_min[1] = loc_u + ind_num;
        }
        else {
            my_min[0] = INFINITY;
            my_min[1] = -1;
        }

        MPI_Allreduce(my_min, glbl_min, 1, MPI_2INT, MPI_MINLOC, comm);

        dist_glbl_u = glbl_min[0];
        glbl_u = glbl_min[1];
        if (glbl_u == -1)
            break;

        if (glbl_u >= ind_num && glbl_u < ind_num+ loc_n) {
            loc_u = glbl_u - ind_num;
            loc_known[loc_u] = 1;
        }

        for (loc_v = 0; loc_v < loc_n; loc_v++) {
            if (!loc_known[loc_v]) {
                new_dist = dist_glbl_u + loc_mat[glbl_u * loc_n + loc_v];
                if (new_dist < loc_dist[loc_v]) {
                    loc_dist[loc_v] = new_dist;
                    loc_pred[loc_v] = glbl_u;
                }
            }
        }
    }
    delete[]loc_known;
}

int Find_min_dist(int loc_dist[], int loc_known[], int loc_n) {
    int loc_u = -1, loc_v;
    int shortest_dist = INFINITY;

    for (loc_v = 0; loc_v < loc_n; loc_v++) {
        if (!loc_known[loc_v]) {
            if (loc_dist[loc_v] < shortest_dist) {
                shortest_dist = loc_dist[loc_v];
                loc_u = loc_v;
            }
        }
    }
    return loc_u;
}

void Print_matrix(int mat[], int rows, int cols) {
    int i, j;

    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++)
            if (mat[i * cols + j] == INFINITY)
                cout << "i ";
            else
                cout <<  mat[i * cols + j]<<" ";
        cout << "\n";
    }

    cout << "\n";
}

void Print_dists(int global_dist[], int n) {
    int v;

    cout <<"  v    dist 0->v\n";
    cout << "----   ---------\n";

    for (v = 1; v < n; v++) {
        if (global_dist[v] == INFINITY) {
            cout << v<< "       inf\n";
        }
        else
            cout << v <<"       " <<global_dist[v]<<"\n";
    }
    cout << "\n";
}


void Print_paths(int global_pred[], int n) {
    int v, w, * path, count, i;

    path = new int[n];

    cout << "  v     Path 0->v\n";
    cout << "----    ---------\n";
    for (v = 1; v < n; v++) {
        cout << v<<":    ";
        count = 0;
        w = v;
        while (w != 0) {
            path[count] = w;
            count++;
            w = global_pred[w];
        }
        cout <<"0 ";
        for (i = count - 1; i >= 0; i--)
            cout << path[i]<<" ";
        cout << "\n";
    }

    delete[]path;
}
