#ifndef __UTILS_H
#define __UTILS_H

#define BIN 0
#define DEC 1
#define DEBUG_MODE 1

typedef enum
{
    VITERBI_HARD,
    VITERBI_SOFT,
    BCJR
} DECODE_METHOD;

typedef struct
{
    int reg_num;
    int *trans_mat;
    int nin;
    int nout;
} PARAMETER;

typedef struct TrellisNode
{
    int data;
    struct TrellisLine *LeftLines;
    struct TrellisLine *RightLines;
} TNode;

typedef struct TrellisLine
{
    int id;
    int input;
    int output;
    struct TrellisNode *BeginNode;
    struct TrellisNode *EndNode;
} TLine;

typedef struct VITERBI_MAT_NODE
{
    float min_cost;
    float min_cost_path;
    int state;
    int active;
} VNODE;

int count_ones(int num);
int *int2array(int num, int len, int mode);
int array2int(int *arr, int len);

PARAMETER *get_essential_params(int A, int B, int C);

#endif