#include "utils.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

int count_ones(int num)
{
    int count = 0;
    while (num)
    {
        count += num & 1;
        num >>= 1;
    }
    return count;
}

// 获得数字的长度，mode0为二进制序列的长度，mode1为十进制的长度
int get_len(int num, int mode)
{
    int len = 0;
    while (num)
    {
        len++;
        switch (mode)
        {
        case 0:
            num >>= 1;
            break;
        case 1:
            num /= 10;
            break;
        }
    }
    return len;
}

int *int2array(int num, int len, int mode)
{
    int *a = (int *)malloc(sizeof(int) * len);
    for (int i = 0; i < len; i++)
    {
        switch (mode)
        {
        case 0:
            a[len - i - 1] = num & 1;
            num >>= 1;
            break;
        case 1:
            a[len - i - 1] = num % 10;
            num /= 10;
            break;
        }
    }
    return a;
}

int array2int(int *arr, int len)
{
    int num = 0;
    for (int i = 0; i < len; i++)
    {
        num <<= 1;
        num += arr[i];
    }
    return num;
}

/**
 * 生成卷积码必要的参数
 */
PARAMETER *get_essential_params(int a, int b, int N)
{

    PARAMETER *params = (PARAMETER *)malloc(sizeof(PARAMETER));
    params->nin = 1;
    params->nout = 2;
    int section = (int)log2(N);

    int *a_declist = int2array(a, get_len(a, DEC), DEC);
    int *b_declist = int2array(b, get_len(b, DEC), DEC);
    int full_section_num = ceil(log10(a)) - 1;
    int single_section_width = floor(log2(a_declist[0])) + 1;
    params->reg_num = full_section_num * section + single_section_width - 1;

    params->trans_mat = (int *)malloc(sizeof(int) * (params->reg_num + 1) * 2);
    int *single_a = int2array(a_declist[0], single_section_width, BIN);
    int *single_b = int2array(b_declist[0], single_section_width, BIN);
    for (int i = 0; i < single_section_width; i++)
    {
        params->trans_mat[i] = single_a[i];
        params->trans_mat[i + params->reg_num + 1] = single_b[i];
    }
    free(single_a);
    free(single_b);

    for (int i = 0; i < full_section_num; i++)
    {
        int *section_a = int2array(a_declist[i + 1], log2(N), BIN);
        int *section_b = int2array(b_declist[i + 1], log2(N), BIN);

        for (int j = 0; j < section; j++)
        {
            params->trans_mat[single_section_width + i * section + j] = section_a[j];
            params->trans_mat[single_section_width + i * section + j + params->reg_num + 1] = section_b[j];
        }
    }

    return params;
}

// int main(){
//     printf("%d", count_ones(5));
// //     int *a = (int*)malloc(sizeof(int)*3);
// //     a[0] = 1;
// //     a[1] = 0;
// //     a[2] = 0;
// //     printf("%d\n", array2int(a, 3));

//     return 0;
// }