#ifndef __UTILS_H
#define __UTILS_H

#define BIN 0
#define DEC 1
#define DEBUG_MODE 1

int* int2array(int num, int len, int mode);
int array2int(int* arr, int len);

typedef struct {
    int reg_num;
    int* trans_mat;
    float code_rate;
} PARAMETER;

PARAMETER* get_essential_params(int A, int B, int C);

typedef struct TrellisNode{
    int data;
    struct TrellisLine* from;
    struct TrellisLine* to;
} TNode;

typedef struct TrellisLine{
    int input;
    int output;
    struct TrellisNode* begin;
    struct TrellisNode* end;
} TLine;


#endif