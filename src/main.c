#include <omp.h>
/***************************************************
Channel Coding Course Work: conolutional codes
This program template has given the message generator, BPSK modulation, AWGN channel model and BPSK demodulation,
you should first determine the encoder structure, then define the message and codeword length, generate the state table, write the convolutional encoder and decoder.

If you have any question, please contact me via e-mail: yanglj39@mail2.sysu.edu.cn
***************************************************/

/** 修改说明：
1. 原代码中 state_num 更改为 reg_num，为寄存器数量；state_num 为寄存器组的状态数量，state_num = 2^reg_num
*/

// #define  _CRT_SECURE_NO_WARNINGS // vs取消scanf报错
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "utils.h"

#define DEBUG_MODE 0
#define DEBUG_STATE_TABLE 0
#define DEBUG_TRELLIS 0
#define DEBUG_VITERBI 1

// 卷积码：0-10，步长1
// turbo码：0-5，步长0.1

static double SNR_START = 0;
static double SNR_FINISH = 10;
static double SNR_STEP = 1;
static int SEQ_NUM = 100;
static int ITERATION_TIMES = 10;

#define message_length 1000000
#define codeword_length (message_length * 2)
DECODE_METHOD decode_method = TURBO;

int CONV_PARAM_1 = 7;
int CONV_PARAM_2 = 5;

void statetable();
void trellis();
void encoder();
void modulation();
void demodulation();
void channel();
void decoder(float SNR_dB);

// functions for turbo
void turbo_encoder(int message_len);
void turbo_decoder(float Snr_dB);
void turbo_statetable();
void shuffle(int *array, int len_of_array);                                                                                                                 // this function is used to shuffle an array to creat the interleaving pattern
void puncture(int *codeword_array, int len_of_message_array);                                                                                               // this function is used to puncture the codeword in turbo encoder
void BCJR_decoder_for_turbo(float SNR_dB, double **priori_prob, double **posteriori_prob, double **extrinsic_prob, int puncture_flag, int interleave_flag); // this function is used to decode the codeword by BCJR algorithm
// because the code rate in turbo may change if puncture_flag = 1, so we need to reconstruct the following functions
void modulation_for_turbo(); // this function is used to modulate the codeword in turbo encoder because the code rate in turbo may change if puncture_flag = 1
void channel_for_turbo();

float code_rate;      // 码率
float code_rate_real; // 真实码率

// channel coefficient
#define pi 3.1415926
#define INF 0xFFFFFFF
double N0, sgm; // 信道噪声

PARAMETER *Parameter;
int reg_num;   // the number of the register of encoder structure
int state_num; // the number of the state of encoder structure
int input_len; // 每次输入的比特数，在758中为1
int input_states;
int *state_table; // state table, the size should be defined yourself
TLine *TLineTable;
TNode *TNodeTable;

int message[message_length], codeword[codeword_length]; // message and codeword，原信息和编码之后的信息
// int codeword2[codeword_length];
int re_codeword[codeword_length]; // the received codeword
int de_message[message_length];   // the decoding message

double tx_symbol[codeword_length][2]; // the transmitted symbols
double rx_symbol[codeword_length][2]; // the received symbols

FILE *fp; // file pointer

// variables for turbo
#define PRINT_TURBO_CODEWORD 0                  // if PRINT_TURBO = 1, the turbo encoder will print the message and the codeword
int turbo_codeword_length = message_length * 3; // the length of codeword in turbo
int turbo_reg_num = 2;                          // the number of the state of register in the turbo encoder structure
int turbo_state_num;                            // the number of the state of turbo encoder structure
int input_state_num = 2;                        // the number of the input of turbo encoder structure
int *random_interleaving_pattern;               // the random interleaving pattern for turbo
int puncture_flag = 0;                          // if puncture_flag = 1, the turbo encoder will puncture the codeword
double **Turbo_Tx_symbol;                       // the transmitted symbols for turbo
double **Turbo_Rx_symbol;                       // the received symbols for turbo
int *turbo_codeword;
int **turbo_state_table;
double sgm_for_turbo;                   // the noise variance for turbo

int main(int argc, char *argv[])
{
    // 译码模式可以传参，前4个参数是START, FINISH, STEP, SEQ_NUM
    if (argc > 1)
    {
        // 参数是浮点数
        SNR_START = atof(argv[1]);
        SNR_FINISH = atof(argv[2]);
        SNR_STEP = atof(argv[3]);
        SEQ_NUM = atoi(argv[4]);

        if (argc > 5)
        {
            if (strcmp(argv[5], "turbo") == 0)
            {
                decode_method = TURBO;
                // 第二个参数可以选择加迭代次数
                if (argc > 6)
                {
                    ITERATION_TIMES = atoi(argv[6]);
                }
            }
            else if (strcmp(argv[5], "viterbi_hard") == 0)
            {
                decode_method = VITERBI_HARD;
                if (argc > 6)
                {
                    CONV_PARAM_1 = atoi(argv[6]);
                    CONV_PARAM_2 = atoi(argv[7]);
                }
            }
            else if (strcmp(argv[5], "viterbi_soft") == 0)
            {
                decode_method = VITERBI_SOFT;
                if (argc > 6)
                {
                    CONV_PARAM_1 = atoi(argv[6]);
                    CONV_PARAM_2 = atoi(argv[7]);
                }
            }
            else if (strcmp(argv[5], "bcjr") == 0)
            {
                decode_method = BCJR;
            }
            else
            {
                printf("Error: unknown decode method\n");
                exit(1);
            }
            printf("Decode method: %s\n", argv[5]);
        }
    }

    code_rate = (float)message_length / (float)codeword_length;
    code_rate_real = (float)(message_length - reg_num) / codeword_length;
    Parameter = get_essential_params(CONV_PARAM_1, CONV_PARAM_2, 8);
    reg_num = Parameter->reg_num; // the number of the register of encoder structure
    state_num = pow(2, reg_num);  // the number of the state of encoder structure
    input_len = 1;                // 每次输入的比特数，在758中为1
    input_states = pow(2, input_len);

    int i;
    float SNR, start, finish;
    long int bit_error, seq, seq_num;
    double BER;
    double progress;

    // generate state table
    if (decode_method == TURBO)
    {
        turbo_statetable();
    }

    statetable();
    trellis();

    // random seed
    srand((int)time(0));

    // //input the SNR and frame number
    // printf("\nEnter start SNR: ");
    // scanf("%f", &start);
    // printf("\nEnter finish SNR: ");
    // scanf("%f", &finish);
    // printf("\nPlease input the number of message: ");
    // scanf("%d", &seq_num);
    start = SNR_START, finish = SNR_FINISH; // 起始和结束的SNR，浮点数，单位为dB
    float SNR_step = SNR_STEP;              // SNR步长
    seq_num = SEQ_NUM;                      // 仿真次数
    fp = fopen("data.txt", "w");
    // 第一行写入message_length，重复次数，turbo译码迭代次数
    fprintf(fp, "%d %d %d %d %d\n", message_length, seq_num, ITERATION_TIMES, CONV_PARAM_1, CONV_PARAM_2);

    for (SNR = start; SNR <= finish; SNR += SNR_step)
    {
        // channel noise
        N0 = (1.0 / code_rate) / pow(10.0, (float)(SNR) / 10.0);
        sgm = sqrt(N0 / 2);

        // turbo channel noise
        if (puncture_flag)
        {
            sgm_for_turbo = sgm; // if puncture_flag = 1, code rate = 1/2
        }
        else
        {
            sgm_for_turbo = sqrt(((1.0 * 3) / pow(10.0, (float)(SNR) / 10.0)) / 2); // if puncture_flag = 0, code rate = 1/3
        }

        bit_error = 0;

        for (seq = 1; seq <= seq_num; seq++)
        {
            // generate binary message randomly
            /****************
            Pay attention that message is appended by 0 whose number is equal to the state of encoder structure.
            ****************/
            // 随机生成信源序列
            for (i = 0; i < message_length - reg_num; i++)
            {
                message[i] = rand() % 2;
            }
            // 信源序列补零
            for (i = message_length - reg_num; i < message_length; i++)
            {
                message[i] = 0;
            }

            // convolutional encoder
            encoder();
            if (decode_method == TURBO)
            {
                turbo_encoder(message_length);
            }

            // BPSK modulation
            modulation();
            if (decode_method == TURBO)
            {
                modulation_for_turbo();
            }

            // AWGN channel
            channel();
            if (decode_method == TURBO)
            {
                channel_for_turbo();
            }

            // BPSK demodulation, it's needed in hard-decision Viterbi decoder
            if (decode_method == VITERBI_HARD)
            {
                demodulation();
            }

            // convolutional decoder
            decoder(SNR);

            // calculate the number of bit error
            for (i = 0; i < message_length; i++)
            {
                if (message[i] != de_message[i])
                    bit_error++;
            }

            progress = (double)(seq * 100) / (double)seq_num;

            // calculate the BER
            BER = (double)bit_error / (double)(message_length * seq);

            // print the intermediate result
            printf("Progress = %5.1f | SNR = %5.1f | Bit Errors = %5.1ld | BER = %E\r", progress, SNR, bit_error, BER);
            if (DEBUG_MODE)
            {
                printf("[DEBUG]:\n");
                printf("信源序列:");
                for (int i = 0; i < message_length; i++)
                {
                    printf("%d ", message[i]);
                }
                printf("\n");
                printf("编码序列:");
                for (int i = 0; i < codeword_length; i++)
                {
                    printf("%d ", codeword[i]);
                }
                printf("\n");
                printf("接收序列:");
                for (int i = 0; i < codeword_length; i++)
                {
                    printf("%d ", re_codeword[i]);
                }
                printf("\n");
                printf("发送符号:");
                for (int i = 0; i < codeword_length; i++)
                {
                    printf("(%f,%f) ", tx_symbol[i][0], tx_symbol[i][1]);
                }
                printf("\n");
                printf("接收符号:");
                for (int i = 0; i < codeword_length; i++)
                {
                    printf("(%f,%f) ", rx_symbol[i][0], rx_symbol[i][1]);
                }
                printf("\n");
                printf("译码序列:");
                for (int i = 0; i < message_length; i++)
                {
                    printf("%d ", de_message[i]);
                }
                printf("\n\n");
            }
        }

        // calculate the BER
        BER = (double)bit_error / (double)(message_length * seq_num);
        // wrire the result into file data.txt
        fprintf(fp, "%f %E\n", SNR, BER);

        printf("Progress = %5.1f | SNR = %5.1f | Bit Errors = %5.1ld | BER = %E\n", progress, SNR, bit_error, BER);
    }
    fclose(fp);
    // system("pause");
    return 0;
}

void statetable()
{
    // 把state_table填满
    // 仿照课件，state_table有state*2行5列
    // 列：IN, Current State, Next State, Out

    state_table = calloc(state_num * input_states * 5, sizeof(int));
    for (int state = 0; state < state_num; state++)
    {
        for (int input = 0; input < input_states; input++)
        {
            int *inputArray = int2array(input, input_len, BIN);
            int *stateArray = int2array(state, reg_num, BIN);

            int line_num = state * input_states + input;
            state_table[5 * line_num] = input;     // IN
            state_table[5 * line_num + 1] = state; // Current State

            int *nextStateArray = (int *)calloc(reg_num, sizeof(int));
            nextStateArray[0] = inputArray[0];
            for (int i = 0; i < reg_num - 1; i++)
            {
                nextStateArray[i + 1] = stateArray[i];
            }

            state_table[5 * line_num + 2] = array2int(nextStateArray, reg_num); // Next State

            int *outArray = (int *)calloc(reg_num, sizeof(int));
            for (int n = 0; n < Parameter->nout; n++)
            {
                outArray[n] = 0;
                outArray[n] ^= inputArray[0] * Parameter->trans_mat[n * reg_num + 0];
                for (int k = 0; k < reg_num; k++)
                {
                    outArray[n] ^= stateArray[k] * Parameter->trans_mat[n * (reg_num + 1) + k + 1];
                }
            }

            state_table[5 * line_num + 3] = array2int(outArray, Parameter->nout); // Out

            state_table[5 * line_num + 4] = line_num + 1; // Line Number

            free(inputArray);
            free(stateArray);
            free(nextStateArray);
            free(outArray);
        }
    }

    // print state_table
    if (DEBUG_MODE && DEBUG_STATE_TABLE)
    {
        printf("[DEBUG]: state_table\n");
        for (int i = 0; i < state_num * input_states; i++)
        {
            for (int j = 0; j < 5; j++)
            {
                printf("%d ", state_table[i * 5 + j]);
            }
            printf("\n");
        }
        printf("\n");
    }
}

void trellis()
{
    TLineTable = (TLine *)calloc(state_num * input_states, sizeof(TLine));
    TNodeTable = (TNode *)calloc(state_num, sizeof(TNode));

    for (int i = 0; i < state_num; i++)
    {
        TNodeTable[i].data = i;
        TNodeTable[i].LeftLines = (TLine *)calloc(input_states, sizeof(TLine));
        TNodeTable[i].RightLines = (TLine *)calloc(input_states, sizeof(TLine));
    }
    int Tlinenum = state_num * input_states;
    for (int i = 0; i < Tlinenum; i++)
    {
        // 每条line就是state table里的一行
        TLineTable[i].id = i;
        TLineTable[i].input = state_table[i * 5];
        TLineTable[i].output = state_table[i * 5 + 3];
        TLineTable[i].BeginNode = &TNodeTable[state_table[i * 5 + 1]];
        TLineTable[i].EndNode = &TNodeTable[state_table[i * 5 + 2]];

        TNodeTable[state_table[i * 5 + 1]].RightLines[TLineTable[i].input] = TLineTable[i];
        // BUG: 相同input会指向同一个状态导致不能用input来区分两条边
        TNodeTable[state_table[i * 5 + 2]].LeftLines[TLineTable[i].BeginNode->data % 2] = TLineTable[i];
    }

    // print trellis
    if (DEBUG_MODE && DEBUG_TRELLIS)
    {
        printf("[DEBUG]: trellis\n");
        for (int i = 0; i < Tlinenum; i++)
        {
            printf("Line %d: in:%d out:%d BeginNode:%d EndNode:%d\n", i, TLineTable[i].input, TLineTable[i].output, TLineTable[i].BeginNode->data, TLineTable[i].EndNode->data);
        }
        printf("\n");
        // for (int i = 0; i < TNodeTable; i++)
        // {
        // 	printf("Node %d: left:%d right%d", i, TNodeTable[i].LeftLines[0].input, TNodeTable[i].LeftLines[0].input);
        // }
    }
}

void encoder()
{
    // convolution encoder, the input is message[] and the output is codeword[]

    // // 先写一个最简单的(7, 5, 8)
    // // state: s0, s1
    // // input: u
    // // output: c1 = u + s0 + s1, c2 = u + s1
    // int s0 = 0, s1 = 0;
    // printf("\n");
    // for (int ii = 0; ii < message_length; ii++){
    // 	codeword2[ii*2] = (message[ii] + s0 + s1) % 2;
    // 	codeword2[ii*2+1] = (message[ii] + s1) % 2;
    // 	// 移位
    // 	s1 = s0;
    // 	s0 = message[ii];
    // 	printf("%d %d ", codeword2[ii*2], codeword2[ii*2+1]);
    // }

    // printf("\n");
    // 寄存器初始化及清零
    int *curStates = (int *)calloc(reg_num + 1, sizeof(int));
    for (int i = 0; i < reg_num + 1; i++)
    {
        curStates[i] = 0;
    }
    // 对 message 第i位编码
    for (int i = 0; i < message_length; i++)
    {
        // 寄存器先移位
        for (int n = reg_num; n > 0; n--)
        {
            curStates[n] = curStates[n - 1];
        }
        curStates[0] = message[i];
        // 再对每一个输出 (c11, c12):
        for (int j = 0; j < Parameter->nout; j++)
        {
            codeword[Parameter->nout * i + j] = 0;
            for (int k = 0; k < reg_num + 1; k++)
            {
                codeword[Parameter->nout * i + j] ^= curStates[k] * Parameter->trans_mat[(reg_num + 1) * j + k];
            }
        }
    }
    free(curStates);
}

void modulation()
{
    // BPSK modulation
    int i;

    // 0 is mapped to (1,0) and 1 is mapped tp (-1,0)
    for (i = 0; i < codeword_length; i++)
    {
        tx_symbol[i][0] = -1 * (2 * codeword[i] - 1);
        tx_symbol[i][1] = 0;
    }
}
void channel()
{
    // AWGN channel
    int i, j;
    double u, r, g;

    for (i = 0; i < codeword_length; i++)
    {
        for (j = 0; j < 2; j++)
        {
            u = (float)rand() / (float)RAND_MAX;
            if (u == 1.0)
                u = 0.999999;
            r = sgm * sqrt(2.0 * log(1.0 / (1.0 - u)));

            u = (float)rand() / (float)RAND_MAX;
            if (u == 1.0)
                u = 0.999999;
            g = (float)r * cos(2 * pi * u);

            rx_symbol[i][j] = tx_symbol[i][j] + g;
        }
    }
}
void demodulation()
{
    int i;
    double d1, d2;
    for (i = 0; i < codeword_length; i++)
    {
        d1 = (rx_symbol[i][0] - 1) * (rx_symbol[i][0] - 1) + rx_symbol[i][1] * rx_symbol[i][1];
        d2 = (rx_symbol[i][0] + 1) * (rx_symbol[i][0] + 1) + rx_symbol[i][1] * rx_symbol[i][1];
        if (d1 < d2)
            re_codeword[i] = 0;
        else
            re_codeword[i] = 1;
    }

    // print re_codeword
    if (DEBUG_MODE && DEBUG_VITERBI)
    {
        printf("[DEBUG]: re_codeword\n");
        for (int i = 0; i < codeword_length; i++)
        {
            printf("%d ", re_codeword[i]);
        }
        printf("\n\n");
    }
}

void decoder_viterbi(int MODE)
{
    int Col = message_length + 1;
    VNODE *VNodeTable = (VNODE *)calloc(state_num * Col, sizeof(VNODE));

    // 初始化
    for (int i = 0; i < state_num; i++)
    {
        for (int j = 0; j < Col; j++)
        {
            int ij = i * Col + j;
            VNodeTable[ij].state = i;
            // 前几列和后几列单独处理
            if (j < reg_num || j >= Col - reg_num)
            {
                VNodeTable[ij].active = 0;
            }
            else
            {
                VNodeTable[ij].active = 1;
            }
            VNodeTable[ij].min_cost = INF;
            VNodeTable[ij].min_cost_path = INF;
        }
    }
    VNodeTable[0].active = 1;       // 第一列只有 0 处是有效节点
    VNodeTable[Col - 1].active = 1; // 最后一列只有 0 处是有效节点
    // 对图进行裁剪，删掉log2(state_num) = reg_num列的不应该出现的节点
    for (int col = 1; col < reg_num; col++)
    {
        for (int row = 0; row < state_num; row++)
        {
            // 前向检查是否有连接到自己的节点
            for (int nout = 0; nout < Parameter->nout; nout++)
            {
                int preid = TNodeTable[row].LeftLines[nout].BeginNode->data;
                if (VNodeTable[preid * Col + col - 1].active)
                {
                    VNodeTable[row * Col + col].active = 1;
                    break;
                }
            }
            // 后向检查是否有连接到自己的节点
            int colbehind = Col - col - 1;
            for (int nout = 0; nout < Parameter->nout; nout++)
            {
                int nextid = TNodeTable[row].RightLines[nout].EndNode->data;
                if (VNodeTable[nextid * Col + colbehind + 1].active)
                {
                    VNodeTable[row * Col + colbehind].active = 1;
                    break;
                }
            }
        }
    }

    // 开始译码
    VNodeTable[0].min_cost = 0;
    VNodeTable[0].min_cost_path = 0;
    // 计算每一列最短路径
    for (int col = 1; col < Col; col++)
    {
        // 计算每条边的差
        double col_mincost = INF; // 这列节点的最小代价
        int decode_output = INF;  // 这列节点的输出
        for (int row = 0; row < state_num; row++)
        {
            int ij = row * Col + col;
            if (!VNodeTable[ij].active)
                continue;
            // 计算从起点到这个节点的所有边中的最短路径
            for (int nout = 0; nout < Parameter->nout; nout++)
            {
                TLine *preLine = &TNodeTable[row].LeftLines[nout];
                int preid = preLine->BeginNode->data;
                int preij = preid * Col + col - 1;
                if (!VNodeTable[preij].active)
                    continue;
                // 计算当前节点的代价
                double cost = 0.0;
                if (MODE == VITERBI_HARD)
                {
                    // 从re_codeword中提取输出
                    int *codeword_section = (int *)calloc(Parameter->nout, sizeof(int));
                    for (int n = 0; n < Parameter->nout; n++)
                    {
                        codeword_section[n] = re_codeword[(col - 1) * Parameter->nout + n];
                    }
                    int codeword_output = array2int(codeword_section, Parameter->nout);
                    cost = VNodeTable[preij].min_cost + count_ones(TLineTable[preLine->id].output ^ codeword_output);
                    free(codeword_section);
                }
                else if (MODE == VITERBI_SOFT)
                {
                    int *line_output = int2array(TLineTable[preLine->id].output, Parameter->nout, BIN);
                    // 2PSK
                    for (int n = 0; n < Parameter->nout; n++)
                    {
                        line_output[n] = -(2 * line_output[n] - 1);
                    }
                    cost = VNodeTable[preij].min_cost;
                    // 计算欧氏距离
                    for (int n = 0; n < Parameter->nout; n++)
                    {
                        cost += pow(rx_symbol[(col - 1) * Parameter->nout + n][0] - (double)line_output[n], 2); // x
                        cost += pow(rx_symbol[(col - 1) * Parameter->nout + n][1] - 0, 2);                      // y
                    }
                    free(line_output);
                }

                if (cost < VNodeTable[ij].min_cost)
                {
                    VNodeTable[ij].min_cost = cost;
                    VNodeTable[ij].min_cost_path = TNodeTable[row].LeftLines[nout].id;
                }
            }
            // if (VNodeTable[ij].min_cost < col_mincost)
            // {
            // 	col_mincost = VNodeTable[ij].min_cost;
            // 	decode_output = VNodeTable[ij].min_cost_path;
            // }
        }
        if (DEBUG_MODE && DEBUG_VITERBI)
        {
            printf("mincost:%.2f output:%d\n", col_mincost, decode_output);
        }
        // de_message[col - 1] = decode_output;
    }

    // 从后往前遍历
    VNODE *node = &VNodeTable[Col - 1]; // 最后一个节点
    for (int col = Col - 1; col > 0; col--)
    {
        for (int inode = 0; inode < Parameter->nout; inode++)
        {
            TLine *line = &TNodeTable[node->state].LeftLines[inode];
            if (line->id == node->min_cost_path)
            {
                de_message[col - 1] = line->input;
                int leftid = line->BeginNode->data;
                node = &VNodeTable[leftid * Col + col - 1];
                break;
            }
        }
    }

    if (DEBUG_MODE && DEBUG_VITERBI)
    {
        printf("[DEBUG]: VNodeTable:(state,active,cost,path)\n");
        for (int i = 0; i < state_num; i++)
        {
            for (int j = 0; j < Col; j++)
            {
                printf("(%.2f %.2f) ",
                       VNodeTable[i * Col + j].min_cost,
                       VNodeTable[i * Col + j].min_cost_path);
            }
            printf("\n");
        }
        printf("\n");
    }
    free(VNodeTable);
}

void decoder_bcjr()
{
    // 输入rx_symbol,输出de_message
    double squared_sigma = N0 / 2;
    double **alpha = calloc(message_length * 4, sizeof(double));
    double **beta = calloc(message_length * 4, sizeof(double));
    for (int i = 0; i < message_length * 4; i++)
    {
        alpha[i] = calloc(4, sizeof(double));
        beta[i] = calloc(4, sizeof(double));
    }

    double gamma_pie_00;
    double gamma_pie_02;
    double gamma_pie_10;
    double gamma_pie_12;
    double gamma_pie_21;
    double gamma_pie_23;
    double gamma_pie_31;
    double gamma_pie_33;

    double p0, p1;
    alpha[0][0] = 1;
    alpha[0][1] = 0;
    alpha[0][2] = 0;
    alpha[0][3] = 0;
    beta[message_length - 1][0] = 1;
    beta[message_length - 1][1] = 0;
    beta[message_length - 1][2] = 0;
    beta[message_length - 1][3] = 0;
    for (int i = 1; i < message_length; i++)
    {
        double ra1 = rx_symbol[2 * i - 2][0];
        double rb1 = rx_symbol[2 * i - 1][0];
        gamma_pie_00 = exp(-1 / N0 * (pow((ra1 - 1), 2) + pow((rb1 - 1), 2))); // v = 00, +1, +1//看卷积码寄存器状态
        gamma_pie_02 = exp(-1 / N0 * (pow((ra1 + 1), 2) + pow((rb1 + 1), 2))); // v = 11, -1, -1
        gamma_pie_10 = exp(-1 / N0 * (pow((ra1 + 1), 2) + pow((rb1 + 1), 2))); // v = 11, -1,-1

        gamma_pie_12 = exp(-1 / N0 * (pow((ra1 - 1), 2) + pow((rb1 - 1), 2))); // v = 00, +1, +1
        gamma_pie_21 = exp(-1 / N0 * (pow((ra1 + 1), 2) + pow((rb1 - 1), 2))); // v = 10, -1, +1

        gamma_pie_23 = exp(-1 / N0 * (pow((ra1 - 1), 2) + pow((rb1 + 1), 2))); // v = 01, +1,-1
        gamma_pie_31 = exp(-1 / N0 * (pow((ra1 - 1), 2) + pow((rb1 + 1), 2))); // v = 01, +1, -1
        gamma_pie_33 = exp(-1 / N0 * (pow((ra1 + 1), 2) + pow((rb1 - 1), 2))); // v = 10, -1, +1
        alpha[i][0] = alpha[i - 1][0] * gamma_pie_00 + alpha[i - 1][1] * gamma_pie_10;
        alpha[i][1] = alpha[i - 1][2] * gamma_pie_21 + alpha[i - 1][3] * gamma_pie_31;
        alpha[i][2] = alpha[i - 1][0] * gamma_pie_02 + alpha[i - 1][1] * gamma_pie_12;
        alpha[i][3] = alpha[i - 1][2] * gamma_pie_23 + alpha[i - 1][3] * gamma_pie_33;
        double alpha_sum = alpha[i][0] + alpha[i][1] + alpha[i][2] + alpha[i][3];
        alpha[i][0] /= alpha_sum;
        alpha[i][1] /= alpha_sum;
        alpha[i][2] /= alpha_sum;
        alpha[i][3] /= alpha_sum;
    }
    for (int i = message_length - 2; i >= 0; i--)
    {
        double ra1 = rx_symbol[2 * i + 2][0];
        double rb1 = rx_symbol[2 * i + 3][0];
        gamma_pie_00 = exp(-1 / N0 * (pow((ra1 - 1), 2) + pow((rb1 - 1), 2))); // v = 00, +1, +1//看卷积码寄存器状态
        gamma_pie_02 = exp(-1 / N0 * (pow((ra1 + 1), 2) + pow((rb1 + 1), 2))); // v = 11, -1, -1
        gamma_pie_10 = exp(-1 / N0 * (pow((ra1 + 1), 2) + pow((rb1 + 1), 2))); // v = 11, -1,-1

        gamma_pie_12 = exp(-1 / N0 * (pow((ra1 - 1), 2) + pow((rb1 - 1), 2))); // v = 00, +1, +1
        gamma_pie_21 = exp(-1 / N0 * (pow((ra1 + 1), 2) + pow((rb1 - 1), 2))); // v = 10, -1, +1

        gamma_pie_23 = exp(-1 / N0 * (pow((ra1 - 1), 2) + pow((rb1 + 1), 2))); // v = 01, +1,-1
        gamma_pie_31 = exp(-1 / N0 * (pow((ra1 - 1), 2) + pow((rb1 + 1), 2))); // v = 01, +1, -1
        gamma_pie_33 = exp(-1 / N0 * (pow((ra1 + 1), 2) + pow((rb1 - 1), 2))); // v = 10, -1, +1
        beta[i][0] = beta[i + 1][0] * gamma_pie_00 + beta[i + 1][2] * gamma_pie_02;
        beta[i][1] = beta[i + 1][0] * gamma_pie_10 + beta[i + 1][2] * gamma_pie_12;
        beta[i][2] = beta[i + 1][1] * gamma_pie_21 + beta[i + 1][3] * gamma_pie_23;
        beta[i][3] = beta[i + 1][1] * gamma_pie_31 + beta[i + 1][3] * gamma_pie_33;
        double beta_sum = beta[i][0] + beta[i][1] + beta[i][2] + beta[i][3]; // 归一化
        beta[i][0] /= beta_sum;
        beta[i][1] /= beta_sum;
        beta[i][2] /= beta_sum;
        beta[i][3] /= beta_sum;
    }
    for (int i = 0; i < message_length; i++)
    {
        double ra1 = rx_symbol[2 * i][0];
        double rb1 = rx_symbol[2 * i + 1][0];
        gamma_pie_00 = exp(-1 / N0 * (pow((ra1 - 1), 2) + pow((rb1 - 1), 2))); // v = 00, +1, +1//看卷积码寄存器状态
        gamma_pie_02 = exp(-1 / N0 * (pow((ra1 + 1), 2) + pow((rb1 + 1), 2))); // v = 11, -1, -1
        gamma_pie_10 = exp(-1 / N0 * (pow((ra1 + 1), 2) + pow((rb1 + 1), 2))); // v = 11, -1,-1

        gamma_pie_12 = exp(-1 / N0 * (pow((ra1 - 1), 2) + pow((rb1 - 1), 2))); // v = 00, +1, +1
        gamma_pie_21 = exp(-1 / N0 * (pow((ra1 + 1), 2) + pow((rb1 - 1), 2))); // v = 10, -1, +1

        gamma_pie_23 = exp(-1 / N0 * (pow((ra1 - 1), 2) + pow((rb1 + 1), 2))); // v = 01, +1,-1
        gamma_pie_31 = exp(-1 / N0 * (pow((ra1 - 1), 2) + pow((rb1 + 1), 2))); // v = 01, +1, -1
        gamma_pie_33 = exp(-1 / N0 * (pow((ra1 + 1), 2) + pow((rb1 - 1), 2))); // v = 10, -1, +1
        p0 = alpha[i][0] * gamma_pie_00 * beta[i][0] + alpha[i][1] * gamma_pie_10 * beta[i][0] + alpha[i][2] * gamma_pie_21 * beta[i][1] + alpha[i][3] * gamma_pie_31 * beta[i][1];
        p1 = alpha[i][0] * gamma_pie_02 * beta[i][2] + alpha[i][1] * gamma_pie_12 * beta[i][2] + alpha[i][2] * gamma_pie_23 * beta[i][3] + alpha[i][3] * gamma_pie_33 * beta[i][3];
        double p_sum = p0 + p1;
        p0 = p0 / p_sum;
        p1 = p1 / p_sum;
        if (p0 > p1)
            de_message[i] = 0;
        else
            de_message[i] = 1;
    }

    for (int i = 0; i < message_length * 4; i++)
    {
        free(alpha[i]);
        free(beta[i]);
    }
    free(alpha);
    free(beta);
}

void decoder(float SNR)
{
    switch (decode_method)
    {
    case VITERBI_HARD:
        decoder_viterbi(VITERBI_HARD);
        break;
    case VITERBI_SOFT:
        decoder_viterbi(VITERBI_SOFT);
        break;
    case BCJR:
        decoder_bcjr();
        break;
    case TURBO:
        turbo_decoder(SNR);
        break;
    }
}

// state table in the slide of the chapter 6
void turbo_statetable()
{
    turbo_state_num = pow(2, turbo_reg_num);
    turbo_state_table = calloc(turbo_state_num * input_state_num, sizeof(int *));
    for (int state = 0; state < turbo_state_num; state++)
    {
        for (int input_state = 0; input_state < input_state_num; input_state++)
        {
            turbo_state_table[state * input_state_num + input_state] = calloc(6, sizeof(int));
            // input_state
            turbo_state_table[state * input_state_num + input_state][0] = input_state;
            // current_state of the first register
            int Current_S1 = (state / 2) % 2;
            turbo_state_table[state * input_state_num + input_state][1] = Current_S1;
            // current_state of the second register
            int Current_S2 = state % 2;
            turbo_state_table[state * input_state_num + input_state][2] = Current_S2;
            // output
            int Output = (input_state + Current_S2) % 2;
            turbo_state_table[state * input_state_num + input_state][3] = Output;
            // next_state of the first register
            turbo_state_table[state * input_state_num + input_state][4] = Output;
            // next_state of the second register
            turbo_state_table[state * input_state_num + input_state][5] = Current_S1;
        }
    }
}

void shuffle(int *array, int length_of_array)
{
    srand(time(NULL));
    for (int i = length_of_array - 1; i > 0; i--)
    {
        int j = rand() % (i + 1);
        int temp = array[i];
        array[i] = array[j];
        array[j] = temp;
    }
    // array[0] = 0;
    // array[1] = 1;
    // array[2] = 2;
    // array[3] = 3;
    // array[4] = 4;
}

void puncture(int *codeword_array, int len_of_message_array)
{
    // puncture the second parity bit when message_index is even, puncture the first parity bit when message_index is odd
    int *punctured_codeword = calloc(len_of_message_array * 2, sizeof(int));
    for (int message_index = 0; message_index < len_of_message_array; message_index++)
    {
        punctured_codeword[message_index * 2] = codeword_array[message_index * 3];
        punctured_codeword[message_index * 2 + 1] = codeword_array[message_index * 3 + (message_index % 2) + 1];
    }
    free(codeword_array);
    codeword_array = punctured_codeword;
}

void turbo_encoder(int message_len)
{
    // message[0] = 1;
    // message[1] = 1;
    // message[2] = 1;
    // message[3] = 0;
    // message[4] = 0;

    int First_Register = 0;
    int Second_Register = 0;
    int table_row_index;
    int codeword_length_for_turbo = message_len * 3;
    turbo_codeword = calloc(codeword_length_for_turbo, sizeof(int)); // 3 is the code rate
    // encode the message
    // first parity bit
    for (int message_index = 0; message_index < message_len; message_index++)
    {
        // bit tail
        if (message_index >= message_len - 2)
        {
            message[message_index] = Second_Register;
        }
        turbo_codeword[message_index * 3] = message[message_index];                          // message bit
        table_row_index = message[message_index] + First_Register * 4 + Second_Register * 2; // the index of the row in the state table , the index is calculated by the input, register1, register2 and their weight
        turbo_codeword[message_index * 3 + 1] = turbo_state_table[table_row_index][3];       // the first parity bit
        Second_Register = First_Register;                                                    // the next state of the second register
        First_Register = turbo_state_table[table_row_index][4];                              // the next state of the first register
    }
    // create the random interleaving pattern
    random_interleaving_pattern = calloc(message_len, sizeof(int));
    for (int i = 0; i < message_len; i++)
    {
        random_interleaving_pattern[i] = i;
    }
    shuffle(random_interleaving_pattern, message_len);
    // interleaving the message
    int *interleaved_message = calloc(message_len, sizeof(int));
    for (int i = 0; i < message_len; i++)
    {
        interleaved_message[i] = message[random_interleaving_pattern[i]];
    }
    First_Register = 0;
    Second_Register = 0;
    for (int message_index = 0; message_index < message_len; message_index++)
    {
        table_row_index = interleaved_message[message_index] + First_Register * 4 + Second_Register * 2; // the index of the row in the state table , the index is calculated by the input, register1, register2 and their weight
        turbo_codeword[message_index * 3 + 2] = turbo_state_table[table_row_index][3];                   // the second parity bit
        Second_Register = First_Register;                                                                // the next state of the second register
        First_Register = turbo_state_table[table_row_index][4];                                          // the next state of the first register
    }
    // puncture the codeword
    if (puncture_flag == 1)
    {
        puncture(turbo_codeword, message_len);
        codeword_length_for_turbo = message_len * 2;
    }
    free(interleaved_message);
    // print the codeword
    if (PRINT_TURBO_CODEWORD)
    {
        printf("[DEBUG]: message\n");
        for (int i = 0; i < message_len; i++)
        {
            printf("%d ", message[i]);
        }
        printf("\n\n");
        printf("[DEBUG]: random_interleaving_pattern\n");
        for (int i = 0; i < message_len; i++)
        {
            printf("%d ", random_interleaving_pattern[i]);
        }
        printf("\n\n");
        printf("[DEBUG]: interleaved_message\n");
        for (int i = 0; i < message_len; i++)
        {
            printf("%d ", interleaved_message[i]);
        }
        printf("\n\n");
        printf("[DEBUG]: turbo_codeword\n");
        for (int i = 0; i < codeword_length_for_turbo; i++)
        {
            printf("%d ", turbo_codeword[i]);
        }
        printf("\n\n");
    }
}

void modulation_for_turbo()
{
    // BPSK modulation
    Turbo_Tx_symbol = calloc(turbo_codeword_length, sizeof(double *));
    for (int i = 0; i < turbo_codeword_length; i++)
    {
        Turbo_Tx_symbol[i] = calloc(2, sizeof(double));
    }
    // 0 is mapped to (1,0) and 1 is mapped tp (-1,0)
    for (int i = 0; i < turbo_codeword_length; i++)
    {
        Turbo_Tx_symbol[i][0] = -1 * (2 * turbo_codeword[i] - 1);
        Turbo_Tx_symbol[i][1] = 0;
    }
    free(turbo_codeword);
}

void channel_for_turbo()
{
    // AWGN channel
    int i, j;
    double u, r, g;
    Turbo_Rx_symbol = calloc(turbo_codeword_length, sizeof(double *));
    for (int i = 0; i < turbo_codeword_length; i++)
    {
        Turbo_Rx_symbol[i] = calloc(2, sizeof(double));
    }

    for (i = 0; i < turbo_codeword_length; i++)
    {
        for (j = 0; j < 2; j++)
        {
            u = (float)rand() / (float)RAND_MAX;
            if (u == 1.0)
                u = 0.999999;
            r = sgm_for_turbo * sqrt(2.0 * log(1.0 / (1.0 - u)));

            u = (float)rand() / (float)RAND_MAX;
            if (u == 1.0)
                u = 0.999999;
            g = (float)r * cos(2 * pi * u);

            Turbo_Rx_symbol[i][j] = Turbo_Tx_symbol[i][j] + g;
        }
    }
    for (int i = 0; i < turbo_codeword_length; i++)
    {
        free(Turbo_Tx_symbol[i]);
    }
    free(Turbo_Tx_symbol);
}

void BCJR_decoder_for_turbo(float SNR_dB, double **priori_prob, double **posteriori_prob, double **extrinsic_prob, int puncture_flag, int interleave_flag)
{
    /**
     * @brief BCJR decoder for turbo code
     * @param Snr the signal noise ratio
     * @param priori_prob the priori probability of the message, priori_prob[0] is the probability of the message is 0, priori_prob[1] is the probability of the message is 1
     * @param posteriori_prob the posteriori probability of the message, posteriori_prob[i][0] is the probability of the message is 0, posteriori_prob[i][1] is the probability of the message is 1
     * @param extrinsic_prob the extrinsic probability of the message, extrinsic_prob[0] is the probability of the message is 0, extrinsic_prob[1] is the probability of the message is 1
     */
    double turbo_code_rate = (double)message_length / turbo_codeword_length;
    // N0 = (1.0 / code_rate) / pow(10.0, (float)(SNR) / 10.0);
    N0 = (1.0 / turbo_code_rate) / pow(10.0, SNR_dB / 10.0);
    sgm = sqrt(N0 / 2);
    double **alpha = calloc(turbo_codeword_length * 4, sizeof(double *));
    double **beta = calloc(turbo_codeword_length * 4, sizeof(double *));
    double *gamma_pie_00 = calloc(message_length, sizeof(double)); // the probability of state transition from 00 to 00
    double *gamma_pie_02 = calloc(message_length, sizeof(double)); // the probability of state transition from 00 to 10
    double *gamma_pie_10 = calloc(message_length, sizeof(double)); // the probability of state transition from 10 to 00
    double *gamma_pie_12 = calloc(message_length, sizeof(double)); // the probability of state transition from 10 to 10
    double *gamma_pie_21 = calloc(message_length, sizeof(double)); // the probability of state transition from 01 to 11
    double *gamma_pie_23 = calloc(message_length, sizeof(double)); // the probability of state transition from 01 to 01
    double *gamma_pie_31 = calloc(message_length, sizeof(double)); // the probability of state transition from 11 to 01
    double *gamma_pie_33 = calloc(message_length, sizeof(double)); // the probability of state transition from 11 to 11

    double **Pch_1 = calloc(message_length, sizeof(double *));            // 信道观察概率  Pch_1[i][0] is the probability of the message bit is 0, Pch_1[i][1] is the probability of the message bit is 1
    double **interleave_Pch_1 = calloc(message_length, sizeof(double *)); // 交织器交织后的信道观察概率？？ Pch_1[i][0] is the probability of the message bit is 0, Pch_1[i][1] is the probability of the message bit is 1
    double **Pch_2 = calloc(message_length, sizeof(double *));            // 校验位的信道观察概率？？？ Pch_2[i][0] is the probability of the parity bit is 0, Pch_2[i][1] is the probability of the parity bit is 1

    #pragma omp parallel for
    for (int i = 0; i < message_length; i++)
    {
        alpha[i] = calloc(4, sizeof(double));
        beta[i] = calloc(4, sizeof(double));
        Pch_1[i] = calloc(2, sizeof(double));
        interleave_Pch_1[i] = calloc(2, sizeof(double));
        Pch_2[i] = calloc(2, sizeof(double));
    }

    double p0, p1; // the probability of the bit is 0 or 1
    // initialize the alpha and beta
    alpha[0][0] = 1; // the probability of begin at the state  00
    alpha[0][1] = 0; // the probability of begin at the state  01
    alpha[0][2] = 0; // the probability of begin at the state  10
    alpha[0][3] = 0; // the probability of begin at the state  11
                     //	beta[message_length - 1][0] = 0; // the probability of ending at the state 00
                     //	beta[message_length - 1][1] = 1; // the probability of ending at the state 01
                     //	beta[message_length - 1][2] = 0; // the probability of ending at the state 10
                     //	beta[message_length - 1][3] = 0; // the probability of ending at the state 11
    if (interleave_flag)
    {
        beta[message_length - 1][0] = 0.25; // the probability of ending at the state 00
        beta[message_length - 1][1] = 0.25; // the probability of ending at the state 01
        beta[message_length - 1][2] = 0.25; // the probability of ending at the state 10
        beta[message_length - 1][3] = 0.25; // the probability of ending at the state 11
    }
    else
    {
        beta[message_length - 1][0] = 1; // the probability of ending at the state 00
        beta[message_length - 1][1] = 0; // the probability of ending at the state 01
        beta[message_length - 1][2] = 0; // the probability of ending at the state 10
        beta[message_length - 1][3] = 0; // the probability of ending at the state 11
    }
    // calculate the Pch for each state
    #pragma omp parallel for
    for (int i = 0; i < message_length; i++)
    {
        // 0 is mapped to (1,0) and 1 is mapped tp (-1,0)
        if (puncture_flag)
        {
            Pch_1[i][0] = exp(-1 / N0 * (pow((Turbo_Rx_symbol[2 * i][0] - 1), 2) + pow(Turbo_Rx_symbol[2 * i][1] - 0, 2))); // output message codeword = 0  which is mapped to (1,0)
            Pch_1[i][1] = exp(-1 / N0 * (pow((Turbo_Rx_symbol[2 * i][0] + 1), 2) + pow(Turbo_Rx_symbol[2 * i][1] - 0, 2))); // output message codeword = 1  which is mapped to (-1,0)
            if (interleave_flag)
            {
                // the odd parity bit is punctured
                if (i % 2 == 0)
                {
                    Pch_2[i][0] = exp(-1 / N0 * (pow((Turbo_Rx_symbol[2 * i + 1][0] - 1), 2) + pow(Turbo_Rx_symbol[2 * i + 1][1] - 0, 2))); // output parity codeword = 0  which is mapped to (1,0)
                    Pch_2[i][1] = exp(-1 / N0 * (pow((Turbo_Rx_symbol[2 * i + 1][0] + 1), 2) + pow(Turbo_Rx_symbol[2 * i + 1][1] - 0, 2))); // output parity codeword = 1  which is mapped to (-1,0)
                }
                else
                {
                    Pch_2[i][0] = 0.5000;
                    Pch_2[i][1] = 0.5000;
                }
            }
            else
            {
                // the even parity bit is punctured
                if (i % 2 == 0)
                {
                    Pch_2[i][0] = 0.5000;
                    Pch_2[i][1] = 0.5000;
                }
                else
                {
                    Pch_2[i][0] = exp(-1 / N0 * (pow((Turbo_Rx_symbol[2 * i + 1][0] - 1), 2) + pow(Turbo_Rx_symbol[2 * i + 1][1] - 0, 2))); // output parity codeword = 0  which is mapped to (1,0)
                    Pch_2[i][1] = exp(-1 / N0 * (pow((Turbo_Rx_symbol[2 * i + 1][0] + 1), 2) + pow(Turbo_Rx_symbol[2 * i + 1][1] - 0, 2))); // output parity codeword = 1  which is mapped to (-1,0)
                }
            }
        }
        else
        {
            if (interleave_flag)
            {
                // message bit
                Pch_1[i][0] = exp(-1 / N0 * (pow((Turbo_Rx_symbol[3 * i][0] - 1), 2) + pow(Turbo_Rx_symbol[3 * i][1] - 0, 2))); // output message codeword = 0  which is mapped to (1,0)
                Pch_1[i][1] = exp(-1 / N0 * (pow((Turbo_Rx_symbol[3 * i][0] + 1), 2) + pow(Turbo_Rx_symbol[3 * i][1] - 0, 2))); // output message codeword = 1  which is mapped to (-1,0)
                // parity bit
                Pch_2[i][0] = exp(-1 / N0 * (pow((Turbo_Rx_symbol[3 * i + 2][0] - 1), 2) + pow(Turbo_Rx_symbol[3 * i + 2][1] - 0, 2))); // output parity codeword = 0  which is mapped to (1,0)
                Pch_2[i][1] = exp(-1 / N0 * (pow((Turbo_Rx_symbol[3 * i + 2][0] + 1), 2) + pow(Turbo_Rx_symbol[3 * i + 2][1] - 0, 2))); // output parity codeword = 1  which is mapped to (-1,0)
            }
            else
            {
                // message bit
                Pch_1[i][0] = exp(-1 / N0 * (pow((Turbo_Rx_symbol[3 * i][0] - 1), 2) + pow(Turbo_Rx_symbol[3 * i][1] - 0, 2))); // output message codeword = 0  which is mapped to (1,0)
                Pch_1[i][1] = exp(-1 / N0 * (pow((Turbo_Rx_symbol[3 * i][0] + 1), 2) + pow(Turbo_Rx_symbol[3 * i][1] - 0, 2))); // output message codeword = 1  which is mapped to (-1,0)
                // parity bit
                Pch_2[i][0] = exp(-1 / N0 * (pow((Turbo_Rx_symbol[3 * i + 1][0] - 1), 2) + pow(Turbo_Rx_symbol[3 * i + 1][1] - 0, 2))); // output parity codeword = 0  which is mapped to (1,0)
                Pch_2[i][1] = exp(-1 / N0 * (pow((Turbo_Rx_symbol[3 * i + 1][0] + 1), 2) + pow(Turbo_Rx_symbol[3 * i + 1][1] - 0, 2))); // output parity codeword = 1  which is mapped to (-1,0)
            }
        }
    }
    // interleaving the probability of the message bit for the second parity bit
    if (interleave_flag)
    {
        #pragma omp parallel for
        for (int i = 0; i < message_length; i++)
        {
            interleave_Pch_1[i][0] = Pch_1[random_interleaving_pattern[i]][0];
            interleave_Pch_1[i][1] = Pch_1[random_interleaving_pattern[i]][1];
        }
        #pragma omp parallel for
        for (int i = 0; i < message_length; i++)
        {
            Pch_1[i][0] = interleave_Pch_1[i][0];
            Pch_1[i][1] = interleave_Pch_1[i][1];
        }
    }

    // calculate the gamma for each state
    #pragma omp parallel for
    for (int i = 0; i < message_length; i++)
    {
        // 0 is mapped to (1,0) and 1 is mapped tp (-1,0)
        gamma_pie_00[i] = priori_prob[i][0] * Pch_1[i][0] * Pch_2[i][0]; // register1 = 0, register2 = 0 -> register1 = 0, register2 = 0  while input = 0,  output = {0, 0}
        gamma_pie_02[i] = priori_prob[i][1] * Pch_1[i][1] * Pch_2[i][1]; // register1 = 0, register2 = 0 -> register1 = 1, register2 = 0  while input = 1,  output = {1, 1}
        gamma_pie_10[i] = priori_prob[i][1] * Pch_1[i][1] * Pch_2[i][0]; // register1 = 0, register2 = 1 -> register1 = 0, register2 = 0  while input = 1,  output = {1, 0}
        gamma_pie_12[i] = priori_prob[i][0] * Pch_1[i][0] * Pch_2[i][1]; // register1 = 0, register2 = 1 -> register1 = 1, register2 = 0  while input = 0,  output = {0, 1}
        gamma_pie_21[i] = priori_prob[i][0] * Pch_1[i][0] * Pch_2[i][0]; // register1 = 1, register2 = 0 -> register1 = 0, register2 = 1  while input = 0,  output = {0, 0}
        gamma_pie_23[i] = priori_prob[i][1] * Pch_1[i][1] * Pch_2[i][1]; // register1 = 1, register2 = 0 -> register1 = 1, register2 = 1  while input = 1,  output = {1, 1}
        gamma_pie_31[i] = priori_prob[i][1] * Pch_1[i][1] * Pch_2[i][0]; // register1 = 1, register2 = 1 -> register1 = 0, register2 = 1  while input = 1,  output = {1, 0}
        gamma_pie_33[i] = priori_prob[i][0] * Pch_1[i][0] * Pch_2[i][1]; // register1 = 1, register2 = 1 -> register1 = 1, register2 = 1  while input = 0,  output = {0, 1}
    }

    // calculate the alpha
    //#pragma omp parallel for
    for (int i = 1; i < message_length; i++)
    {
        alpha[i][0] = alpha[i - 1][0] * gamma_pie_00[i - 1] + alpha[i - 1][1] * gamma_pie_10[i - 1];
        alpha[i][1] = alpha[i - 1][2] * gamma_pie_21[i - 1] + alpha[i - 1][3] * gamma_pie_31[i - 1];
        alpha[i][2] = alpha[i - 1][0] * gamma_pie_02[i - 1] + alpha[i - 1][1] * gamma_pie_12[i - 1];
        alpha[i][3] = alpha[i - 1][2] * gamma_pie_23[i - 1] + alpha[i - 1][3] * gamma_pie_33[i - 1];
        // normalize
        double alpha_sum = alpha[i][0] + alpha[i][1] + alpha[i][2] + alpha[i][3];
        alpha[i][0] /= alpha_sum;
        alpha[i][1] /= alpha_sum;
        alpha[i][2] /= alpha_sum;
        alpha[i][3] /= alpha_sum;
    }

    // calculate the beta
    //#pragma omp parallel for
    for (int i = message_length - 2; i >= 0; i--)
    {
        beta[i][0] = beta[i + 1][0] * gamma_pie_00[i + 1] + beta[i + 1][2] * gamma_pie_02[i + 1];
        beta[i][1] = beta[i + 1][0] * gamma_pie_10[i + 1] + beta[i + 1][2] * gamma_pie_12[i + 1];
        beta[i][2] = beta[i + 1][1] * gamma_pie_21[i + 1] + beta[i + 1][3] * gamma_pie_23[i + 1];
        beta[i][3] = beta[i + 1][1] * gamma_pie_31[i + 1] + beta[i + 1][3] * gamma_pie_33[i + 1];

        //		beta[i][1] = beta[i + 1][0] * gamma_pie_10[i] + beta[i + 1][2] * gamma_pie_12[i];
        //		beta[i][2] = beta[i + 1][1] * gamma_pie_21[i] + beta[i + 1][3] * gamma_pie_23[i];
        //		beta[i][3] = beta[i + 1][1] * gamma_pie_31[i] + beta[i + 1][3] * gamma_pie_33[i];
        // normalize
        double beta_sum = beta[i][0] + beta[i][1] + beta[i][2] + beta[i][3];
        beta[i][0] /= beta_sum;
        beta[i][1] /= beta_sum;
        beta[i][2] /= beta_sum;
        beta[i][3] /= beta_sum;
    }

    // calculate the posteriori probability
    #pragma omp parallel for
    for (int i = 0; i < message_length; i++)
    {
        // posteriori_prob = sigma {alpha * gamma * beta}
        posteriori_prob[i][0] = alpha[i][0] * gamma_pie_00[i] * beta[i][0] + alpha[i][1] * gamma_pie_12[i] * beta[i][2] + alpha[i][2] * gamma_pie_21[i] * beta[i][1] + alpha[i][3] * gamma_pie_33[i] * beta[i][3];
        posteriori_prob[i][1] = alpha[i][0] * gamma_pie_02[i] * beta[i][2] + alpha[i][1] * gamma_pie_10[i] * beta[i][0] + alpha[i][2] * gamma_pie_23[i] * beta[i][3] + alpha[i][3] * gamma_pie_31[i] * beta[i][1];
        // normalize
        double posteriori_prob_sum = posteriori_prob[i][0] + posteriori_prob[i][1];
        posteriori_prob[i][0] /= posteriori_prob_sum;
        posteriori_prob[i][1] /= posteriori_prob_sum;
    }

    // calculate the extrinsic probability
    #pragma omp parallel for
    for (int i = 0; i < message_length; i++)
    {
        // extrinsic_prob = posteriori_prob / priori_prob
        extrinsic_prob[i][0] = posteriori_prob[i][0] / priori_prob[i][0];
        extrinsic_prob[i][1] = posteriori_prob[i][1] / priori_prob[i][1];
        // normalize
        double extrinsic_prob_sum = extrinsic_prob[i][0] + extrinsic_prob[i][1];
        extrinsic_prob[i][0] /= extrinsic_prob_sum;
        extrinsic_prob[i][1] /= extrinsic_prob_sum;
    }

    // free the memory
    #pragma omp parallel for
    for (int i = 0; i < message_length; i++)
    {
        free(alpha[i]);
        free(beta[i]);
        free(Pch_1[i]);
        free(interleave_Pch_1[i]);
        free(Pch_2[i]);
    }
    free(alpha);
    free(beta);
    free(gamma_pie_00);
    free(gamma_pie_02);
    free(gamma_pie_10);
    free(gamma_pie_12);
    free(gamma_pie_21);
    free(gamma_pie_23);
    free(gamma_pie_31);
    free(gamma_pie_33);
    free(Pch_1);
    free(interleave_Pch_1);
    free(Pch_2);
}

void turbo_decoder(float Snr_dB)
{
    double **Priori_prob = calloc(message_length, sizeof(double *));
    double **Posteriori_prob = calloc(message_length, sizeof(double *));
    double **Extrinsic_prob = calloc(message_length, sizeof(double *));
    double **Posteriori_prob_final = calloc(message_length, sizeof(double *));
    #pragma omp parallel for
    for (int i = 0; i < message_length; i++)
    {
        Priori_prob[i] = calloc(2, sizeof(double));
        Posteriori_prob[i] = calloc(2, sizeof(double));
        Extrinsic_prob[i] = calloc(2, sizeof(double));
        Posteriori_prob_final[i] = calloc(2, sizeof(double));
    }
    // initialize the priori probability
    #pragma omp parallel for
    for (int i = 0; i < message_length; i++)
    {
        Priori_prob[i][0] = 0.5;
        Priori_prob[i][1] = 0.5;
    }
    // turbo decoding
    // #pragma omp parallel for
    for (int i = 0; i < ITERATION_TIMES; i++)
    {
        BCJR_decoder_for_turbo(Snr_dB, Priori_prob, Posteriori_prob, Extrinsic_prob, puncture_flag, 0);
        // interleaving the probabilities for the BCJR(2)
        #pragma omp parallel for
        for (int j = 0; j < message_length; j++)
        {
            Priori_prob[j][0] = Extrinsic_prob[random_interleaving_pattern[j]][0];
            Priori_prob[j][1] = Extrinsic_prob[random_interleaving_pattern[j]][1];
        }

        if (PRINT_TURBO_CODEWORD)
        {
            printf("[DEBUG]: Priori_prob\n");
            for (int i = 0; i < message_length; i++)
            {
                printf("%lf %lf\n", Priori_prob[i][0], Priori_prob[i][1]);
            }
            printf("\n");
            printf("[DEBUG]: Posteriori_prob\n");
            for (int i = 0; i < message_length; i++)
            {
                printf("%lf %lf\n", Posteriori_prob[i][0], Posteriori_prob[i][1]);
            }
            printf("\n");
            printf("[DEBUG]: Extrinsic_prob\n");
            for (int i = 0; i < message_length; i++)
            {
                printf("%lf %lf\n", Extrinsic_prob[i][0], Extrinsic_prob[i][1]);
            }
            printf("\n");
        }

        BCJR_decoder_for_turbo(Snr_dB, Priori_prob, Posteriori_prob, Extrinsic_prob, puncture_flag, 1);
        // inverse interleaving the probabilities for the BCJR(1)
        #pragma omp parallel for
        for (int j = 0; j < message_length; j++)
        {
            Priori_prob[random_interleaving_pattern[j]][0] = Extrinsic_prob[j][0];
            Priori_prob[random_interleaving_pattern[j]][1] = Extrinsic_prob[j][1];
        }
    }
    // inverse the posteriori to get the final posteriori probability
    #pragma omp parallel for
    for (int i = 0; i < message_length; i++)
    {
        Posteriori_prob_final[random_interleaving_pattern[i]][0] = Posteriori_prob[i][0];
        Posteriori_prob_final[random_interleaving_pattern[i]][1] = Posteriori_prob[i][1];
    }
    // compare the posteriori probability to get the decoded message
    #pragma omp parallel for
    for (int i = 0; i < message_length; i++)
    {
        if (Posteriori_prob_final[i][0] > Posteriori_prob_final[i][1])
        {
            de_message[i] = 0;
        }
        else
        {
            de_message[i] = 1;
        }
    }
    // free the memory
    #pragma omp parallel for
    for (int i = 0; i < message_length; i++)
    {
        free(Priori_prob[i]);
        free(Posteriori_prob[i]);
        free(Extrinsic_prob[i]);
        free(Posteriori_prob_final[i]);
    }
    free(Priori_prob);
    free(Posteriori_prob);
    free(Extrinsic_prob);
    free(Posteriori_prob_final);
    for (int i = 0; i < turbo_codeword_length; i++)
    {
        free(Turbo_Rx_symbol[i]);
    }
    free(Turbo_Rx_symbol);
    free(random_interleaving_pattern);
}