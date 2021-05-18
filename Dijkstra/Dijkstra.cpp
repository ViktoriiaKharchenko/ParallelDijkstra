// Dijkstra.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
#include <time.h>
#include <iostream>
#include <chrono>
#define INFINITY 10000000
#include<iomanip>

using namespace std;
void Read_matrix(int mat[], int n);
void Dijkstra_Init(int mat[], int pred[], int dist[], int known[], int n);
void RandomDataInitialization(int* pMatrix, int n);
void Dijkstra(int mat[], int dist[], int pred[], int n);
int Find_min_dist(int dist[], int known[], int n);
void Print_matrix(int mat[], int rows, int cols);
void Print_dists(int dist[], int n);
void Print_paths(int pred[], int n);

int main()
{
    int *dist, *pred, *mat;
    int n;
    time_t  start, finish;
    double duration;
    cin >> n;  
    mat = new int[n * n];
    dist = new int[n];
    pred = new int[n];
    //Read_matrix(mat,n);

    RandomDataInitialization(mat, n);
    start = clock();

    Dijkstra(mat, dist, pred, n);
    finish = clock();
    duration = (finish - start) / double(CLOCKS_PER_SEC);
    cout <<"\n Time of execution: " <<duration;
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

void Dijkstra_Init(int mat[], int pred[], int dist[], int known[], int n) {
    for (int i= 0; i < n; i++) {
        known[i] = 0;
        dist[i] = INFINITY;
        pred[i] = 0;
    }
}

void Dijkstra(int mat[], int dist[], int pred[], int n) {

    int  v, u,  new_dist;
    int* known= new int[n];
    int my_min;

    Dijkstra_Init(mat, pred, dist, known, n);
    dist[0] = 0;

    for (v = 0; v < n-1; v++) {
        u = Find_min_dist(dist, known, n);
        if (u == -1) break;
        known[u] = 1;
        if (u != -1) {
            my_min = dist[u];
        }
        else {
            my_min = INFINITY;
        }

        for (int i = 0; i < n; i++) {
            if (!known[i]) {
                new_dist = my_min + mat[u * n + i];
                if (new_dist < dist[i]) {
                    dist[i] = new_dist;
                    pred[i] = u;
                }
            }
        }
    }
    delete[]known;
}

int Find_min_dist(int dist[], int known[], int n) {
    int shortest_dist = INFINITY;
    int u = -1;
    for (int i = 0; i < n; i++) {
        if (!known[i]) {
            if (dist[i] < shortest_dist) {
                shortest_dist = dist[i];
                u = i;
            }
        }
    }
    return u;
}
void Print_matrix(int mat[], int rows, int cols) {
    int i, j;

    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++)
            if (mat[i * cols + j] == INFINITY)
                cout <<"i ";
            else
                cout << mat[i * cols + j]<<" ";
        cout << "\n";
    }

    cout << "\n";
}


void Print_dists(int dist[], int n) {
    int v;

    cout <<"  v    dist 0->v\n";
    cout <<"----   ---------\n";

    for (v = 1; v < n; v++) {
        if (dist[v] == INFINITY) {
            cout <<v << "       inf\n";
        }
        else
            cout << v <<"       " << dist[v]<<"\n";
    }
     cout << "\n";
}


void Print_paths(int pred[], int n) {
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
            w = pred[w];
        }
        cout << "0 ";
        for (i = count - 1; i >= 0; i--)
            cout << path[i] <<" ";
        cout << "\n";
    }
    delete[]path;
}