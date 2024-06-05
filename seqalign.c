#include <stdio.h>
#include <string.h>
#include <stdlib.h>

void printAnno(){

 printf("*---------------------------------------------------------*\n\n\
          @Func: This program is to implement a simple sequence alignment;\n\
          @Course:         Bioinfo methods;\n\
          @Teacher:        Ruan Jue;\n\
          @Date:           2022/9/5;\n\
          @Code Modifier:  Chen S\n\n\
          @Scoring rules:      Match = 1, Mismatch = -3, Extend = -2;\n\
         *---------------------------------------------------------*\n\n");
}


// Initialize the alignment matrix
void init_align_matrix(int ql, int tl, int **mtx, int E){

    int x,y;
    mtx[0][0] = 0;
    for(x = 1; x <= ql; x++) mtx[0][x] = mtx[0][x-1] + E;
    for(y = 1; y <= tl; y++) mtx[y][0] = mtx[y-1][0] + E;
}

//Get max value from three numbers 
int n_max(int x, int y){

    return (x > y) ? x : y;
}

int num_max(int a, int b, int c){
    int max;
    max = n_max(a, n_max(b,c));
    return max;
}


//generate matrix
/*
int **generate_matrix(int rows,int columns)
{
    int **numbers=new int*[rows];
    for(int i=0;i<rows;i++){
        numbers[i]=new int[columns];
        for(int j=0;j<columns;j++)
                numbers[i][j]=i*columns+j;
    }
    return numbers;
}
*/

int **generate_matrix(int rows,int columns)
{
    int **numbers= (int**)malloc(sizeof (int *) * rows);
    for(int i=0;i<rows;i++){
        numbers[i]= (int *)malloc(sizeof (int) * columns);
        for(int j=0;j<columns;j++)
                numbers[i][j]=i*columns+j;
    }
    return numbers;
}

void calc_align_matrix(char *qseq, int ql, char *tseq, int tl, int **mtx, int M, int X, int E){

    int x, y;
    for(y=1; y<=tl; y++){
        for(x=1; x<=ql; x++){
            mtx[x][y] = num_max( mtx[x-1][y-1] + (qseq[x-1] == tseq[y-1]? M : X), mtx[x-1][y] + E, mtx[x][y-1] + E );
            if(x == ql && y== tl) printf("The alignment score is : %d\n", mtx[ql][tl]);
        }
    }
}


void trac_align_matrix(char *qseq, int ql, char *tseq, int tl, int **mtx, int M, int X, int E){

    char qstr[201], tstr[201], cstr[201];
    int x, y, cl;
    x = ql, y = tl; cl = 0;
    while(x > 0 && y > 0 ){
        if(mtx[x][y] == mtx[x-1][y-1] + (qseq[x-1] == tseq[y-1]? M : X)){
            qstr[cl] = qseq[x-1];
            tstr[cl] = tseq[y-1];
            cstr[cl++] = (qseq[--x] == tseq [--y]? '|' : '*' );
        }
        else if(mtx[x][y] == mtx[x-1][y] + E){
            qstr[cl] = qseq[--x];
            tstr[cl] = '-';
            cstr[cl++] = '-';
        }
        else{
            tstr[cl] = tseq[--y];
            qstr[cl] = '-';
            cstr[cl++] = '-';
        }
    }

    for(x=cl-1;x>=0;x--) putchar(qstr[x]); putchar('\n');
    for(x=cl-1;x>=0;x--) putchar(cstr[x]); putchar('\n');
    for(x=cl-1;x>=0;x--) putchar(tstr[x]); putchar('\n');
}


int main(int argc, char **argv){

    printAnno();

    char qseq[101], tseq[101];   //Define two array
    int **mtx = NULL;
    int M, X, E;
    int x, y ,ql, tl;
    //match = 1, mismatch = -3, extend = -2;
    M = 1; X = -3; E = -2;
    printf("Please input two seqs in lines, each less than 100 bp:\n");
    scanf("%s", qseq); ql = strlen(qseq);
    scanf("%s", tseq); tl = strlen(tseq);
    printf("Sequences load completed, begining sequence alignment!\n......\n");
    
    mtx = generate_matrix(101,101);

    init_align_matrix(ql, tl, mtx, E);
    calc_align_matrix(qseq, ql, tseq, tl, mtx, M, X, E);
    trac_align_matrix(qseq, ql, tseq, tl, mtx, M, X, E);

    //free memory
/*    for(int i = 0;i <= 101; i++)
         delete [] mtx[i];
    delete mtx; 
*/
    free(mtx);
    return 0;
}
