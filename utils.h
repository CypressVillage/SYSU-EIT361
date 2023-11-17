#ifndef __UTILS_H
#define __UTILS_H

#define BIN 0
#define DEC 1
#define DEBUG_MODE 1

typedef enum
{
    VICTERBI_HARD,
    VICTERBI_SOFT,
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
    struct TrellisLine *from;
    struct TrellisLine *to;
} TNode;

typedef struct TrellisLine
{
    // int id;
    int input;
    int output;
    struct TrellisNode *begin;
    struct TrellisNode *end;
} TLine;

typedef struct VICTERBI_MAT_NODE
{
    float min_cost;
    float min_cost_path;
    int state;
    int active;
} VNODE;

int *int2array(int num, int len, int mode);
int array2int(int *arr, int len);

PARAMETER *get_essential_params(int A, int B, int C);

#endif