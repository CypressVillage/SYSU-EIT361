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
#include <time.h>
#include <math.h>
#include "utils.h"

void statetable();
void trellis();
void encoder();
void modulation();
void demodulation();
void channel();
void decoder();

#define message_length 1000				   // the length of message
#define codeword_length message_length * 2 // the length of codeword
float code_rate;						   // 码率

// channel coefficient
#define pi 3.1415926
#define INF 0xFFFFFFF
// #define INF 9
double N0, sgm; // 信道噪声

PARAMETER *Parameter;
DECODE_METHOD decode_method = BCJR;
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
int de_message[message_length];	  // the decoding message

double tx_symbol[codeword_length][2]; // the transmitted symbols
double rx_symbol[codeword_length][2]; // the received symbols

FILE *fp; // file pointer

int main()
{
	code_rate = (float)message_length / (float)codeword_length;
	Parameter = get_essential_params(7, 5, 8);
	reg_num = Parameter->reg_num; // the number of the register of encoder structure
	state_num = pow(2, reg_num);  // the number of the state of encoder structure
	input_len = 1;				  // 每次输入的比特数，在758中为1
	input_states = pow(2, input_len);

	int i;
	float SNR, start, finish;
	long int bit_error, seq, seq_num;
	double BER;
	double progress;

	// generate state table
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
	start = -5, finish = 5; // 起始和结束的SNR，浮点数，单位为dB
	float SNR_step = 1;		// SNR步长
	seq_num = 2;			// 仿真次数
	fp = fopen("data.txt", "w");

	for (SNR = start; SNR <= finish; SNR += SNR_step)
	{
		// channel noise
		N0 = (1.0 / code_rate) / pow(10.0, (float)(SNR) / 10.0); // TODO：为什么是1/code_rate
		sgm = sqrt(N0 / 2);

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

			// BPSK modulation
			modulation();

			// AWGN channel
			channel(); // TODO: 记得改回来

			// BPSK demodulation, it's needed in hard-decision Viterbi decoder
			if (decode_method == VITERBI_HARD)
			{
				demodulation();
			}

			// convolutional decoder
			decoder();

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
			printf("Progress=%2.1f, SNR=%2.1f, Bit Errors=%2.1ld, BER=%E\n", progress, SNR, bit_error, BER);
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
				printf("\n\n");
			}
		}

		// calculate the BER
		BER = (double)bit_error / (double)(message_length * seq_num);
		// wrire the result into file data.txt
		fprintf(fp, "%f %E\n", SNR, BER);

		// printf("Progress=%2.1f, SNR=%2.1f, Bit Errors=%2.1d, BER=%E\n", progress, SNR, bit_error, BER);
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
			state_table[5 * line_num] = input;	   // IN
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
	if (DEBUG_MODE)
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
	if (DEBUG_MODE)
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

	// 先写一个最简单的(7, 5, 8)
	// state: s0, s1
	// input: u
	// output: c1 = u + s0 + s1, c2 = u + s1
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
			// printf("%d ", codeword[Parameter->nout * i + j]);
		}
	}
	// printf("\n");
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
	if (DEBUG_MODE)
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
	// codeword[0] = 0;
	// codeword[1] = 0;
	// codeword[2] = 1;
	// codeword[3] = 1;
	// codeword[4] = 1;
	// codeword[5] = 1;
	// codeword[6] = 0;
	// codeword[7] = 1;
	// codeword[8] = 0;
	// codeword[9] = 0;
	// codeword[10] = 1;
	// codeword[11] = 0;
	// codeword[12] = 1;
	// codeword[13] = 1;
	// codeword[14] = 0;
	// codeword[15] = 0;
	// codeword[16] = 1;
	// codeword[17] = 0;
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
	VNodeTable[0].active = 1;		// 第一列只有 0 处是有效节点
	VNodeTable[Col - 1].active = 1; // 最后一列只有 0 处是有效节点
	// 对图进行裁剪，删掉log2(state_num) = reg_num列的不应该出现的节点
	for (int col = 1; col < reg_num; col++)
	{
		for (int row = 0; row < state_num; row++)
		{
			int is_parent_active = 0;
			// 前向检查是否有连接到自己的节点
			for (int nout = 0; nout < Parameter->nout; nout++)
			{
				int preid = TNodeTable[row].LeftLines[nout].BeginNode->data;
				if (VNodeTable[preid * Col + col - 1].active)
				{
					is_parent_active = 1;
				}
			}
			if (is_parent_active)
			{
				VNodeTable[row * Col + col].active = 1;
			}
			// 后向检查是否有连接到自己的节点
			int colbehind = Col - col - 1;
			int is_child_active = 0;
			for (int nout = 0; nout < Parameter->nout; nout++)
			{
				int nextid = TNodeTable[row].RightLines[nout].EndNode->data;
				if (VNodeTable[nextid * Col + colbehind + 1].active)
				{
					is_child_active = 1;
				}
			}
			if (is_child_active)
			{
				VNodeTable[row * Col + colbehind].active = 1;
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
			// 计算到这个节点的所有边中的最短路径
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
					// 从codeword中提取输出
					int *codeword_section = (int *)calloc(Parameter->nout, sizeof(int));
					for (int n = 0; n < Parameter->nout; n++)
					{
						codeword_section[n] = re_codeword[(col - 1) * Parameter->nout + n];
					}
					int codeword_output = array2int(codeword_section, Parameter->nout);
					cost = VNodeTable[preij].min_cost + count_ones(TLineTable[preLine->id].output ^ codeword_output);
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
						cost += pow(rx_symbol[(col - 1) * Parameter->nout + n][1] - 0, 2);						// y
					}
				}

				if (cost < VNodeTable[ij].min_cost)
				{
					VNodeTable[ij].min_cost = cost;
					VNodeTable[ij].min_cost_path = TNodeTable[row].LeftLines[nout].input;
				}
			}
			if (VNodeTable[ij].min_cost < col_mincost)
			{
				col_mincost = VNodeTable[ij].min_cost;
				decode_output = VNodeTable[ij].min_cost_path;
			}
		}
		if (DEBUG_MODE)
		{
			printf("mincost:%.2f output:%d\n", col_mincost, decode_output);
		}
		de_message[col - 1] = decode_output;
	}

	// 后向遍历，找到最短路径
	// int mincost_path = INF;
	// double mincost = INF;
	// for (int col = Col-1; col > 0; col++)
	// {
	// 	for (int row = 0; row < state_num; row++)
	// 	{
	// 		int ij = row * Col + Col - 1;
	// 		if (!VNodeTable[ij].active)
	// 			continue;
	// 		if (VNodeTable[ij].min_cost < mincost)
	// 		{
	// 			mincost = VNodeTable[ij].min_cost;
	// 			mincost_path = VNodeTable[ij].min_cost_path;
	// 		}
	// 	}
	// 	de_message[Col - 1] = mincost_path;
	// }

	if (DEBUG_MODE)
	{
		printf("[DEBUG]: VNodeTable:(state,active,cost,path)\n");
		for (int i = 0; i < state_num; i++)
		{
			for (int j = 0; j < Col; j++)
			{
				printf("(%.2f %.2f) ",
					   //    VNodeTable[i * message_length + j].state,
					   //    VNodeTable[i * message_length + j].active,
					   VNodeTable[i * Col + j].min_cost,
					   VNodeTable[i * Col + j].min_cost_path + 1);
			}
			printf("\n");
		}
		printf("\n");
	}
}

void decoder_bcjr()
{
	// BCJR decoder
	// 输入rx_symbol,输出de_message
	double squared_sigma = N0 / 2;
	double alpha[message_length][4];
	double beta[message_length][4];
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
		// double ra2 = rx_symbol[2 * i-2][1];
		double rb1 = rx_symbol[2 * i - 1][0];
		// double rb2 = rx_symbol[2 * i - 1][1];
		/*gamma_pie_00 = exp(-1 / (2 * squared_sigma) * (pow((ra1 - 1), 2) + pow(ra2, 2) + pow((rb1 - 1), 2) + pow(rb2, 2)));   //v = 00, +1, +1
		gamma_pie_02 = exp(-1 / (2 * squared_sigma) * (pow((ra1 + 1), 2) + pow(ra2, 2) + pow((rb1 + 1), 2) + pow(rb2, 2)));  // v = 11, -1, -1
		gamma_pie_10 = exp(-1 / (2 * squared_sigma) * (pow((ra1 + 1), 2) + pow(ra2, 2) + pow((rb1 - 1), 2) + pow(rb2, 2)));  // v = 10, -1, +1
		gamma_pie_12 = exp(-1 / (2 * squared_sigma) * (pow((ra1 - 1), 2) + pow(ra2, 2) + pow((rb1 + 1), 2) + pow(rb2, 2)));  // v = 01, +1, -1
		gamma_pie_21 = exp(-1 / (2 * squared_sigma) * (pow((ra1 - 1), 2) + pow(ra2, 2) + pow((rb1 - 1), 2) + pow(rb2, 2)));  // v = 00, +1, +1
		gamma_pie_23 = exp(-1 / (2 * squared_sigma) * (pow((ra1 + 1), 2) + pow(ra2, 2) + pow((rb1 + 1), 2) + pow(rb2, 2)));  // v = 11, -1,-1
		gamma_pie_31 = exp(-1 / (2 * squared_sigma) * (pow((ra1 + 1), 2) + pow(ra2, 2) + pow((rb1 - 1), 2) + pow(rb2, 2)));   // v = 10, -1, +1
		gamma_pie_33 = exp(-1 / (2 * squared_sigma) * (pow((ra1- 1), 2) + pow(ra2, 2) + pow((rb1 + 1), 2) + pow(rb2, 2)));   // v = 01, +1, -1*/
		gamma_pie_00 = exp(-1 / (2 * squared_sigma) * (pow((ra1 - 1), 2) + pow((rb1 - 1), 2))); // v = 00, +1, +1//看卷积码寄存器状态
		gamma_pie_02 = exp(-1 / (2 * squared_sigma) * (pow((ra1 + 1), 2) + pow((rb1 + 1), 2))); // v = 11, -1, -1
		gamma_pie_10 = exp(-1 / (2 * squared_sigma) * (pow((ra1 + 1), 2) + pow((rb1 + 1), 2))); // v = 11, -1,-1

		gamma_pie_12 = exp(-1 / (2 * squared_sigma) * (pow((ra1 - 1), 2) + pow((rb1 - 1), 2))); // v = 00, +1, +1
		gamma_pie_21 = exp(-1 / (2 * squared_sigma) * (pow((ra1 + 1), 2) + pow((rb1 - 1), 2))); // v = 10, -1, +1

		gamma_pie_23 = exp(-1 / (2 * squared_sigma) * (pow((ra1 - 1), 2) + pow((rb1 + 1), 2))); // v = 01, +1,-1
		gamma_pie_31 = exp(-1 / (2 * squared_sigma) * (pow((ra1 - 1), 2) + pow((rb1 + 1), 2))); // v = 01, +1, -1
		gamma_pie_33 = exp(-1 / (2 * squared_sigma) * (pow((ra1 + 1), 2) + pow((rb1 - 1), 2))); // v = 10, -1, +1
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
		// double ra2 = rx_symbol[2 * i +2][1];
		double rb1 = rx_symbol[2 * i + 3][0];
		// double rb2 = rx_symbol[2 * i +3][1];
		gamma_pie_00 = exp(-1 / (2 * squared_sigma) * (pow((ra1 - 1), 2) + pow((rb1 - 1), 2))); // v = 00, +1, +1//看卷积码寄存器状态
		gamma_pie_02 = exp(-1 / (2 * squared_sigma) * (pow((ra1 + 1), 2) + pow((rb1 + 1), 2))); // v = 11, -1, -1
		gamma_pie_10 = exp(-1 / (2 * squared_sigma) * (pow((ra1 + 1), 2) + pow((rb1 + 1), 2))); // v = 11, -1,-1

		gamma_pie_12 = exp(-1 / (2 * squared_sigma) * (pow((ra1 - 1), 2) + pow((rb1 - 1), 2))); // v = 00, +1, +1
		gamma_pie_21 = exp(-1 / (2 * squared_sigma) * (pow((ra1 + 1), 2) + pow((rb1 - 1), 2))); // v = 10, -1, +1

		gamma_pie_23 = exp(-1 / (2 * squared_sigma) * (pow((ra1 - 1), 2) + pow((rb1 + 1), 2))); // v = 01, +1,-1
		gamma_pie_31 = exp(-1 / (2 * squared_sigma) * (pow((ra1 - 1), 2) + pow((rb1 + 1), 2))); // v = 01, +1, -1
		gamma_pie_33 = exp(-1 / (2 * squared_sigma) * (pow((ra1 + 1), 2) + pow((rb1 - 1), 2))); // v = 10, -1, +1
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
		// double ra2 = rx_symbol[2 * i][1];
		double rb1 = rx_symbol[2 * i + 1][0];
		// double rb2 = rx_symbol[2 * i + 1][1];
		gamma_pie_00 = exp(-1 / (2 * squared_sigma) * (pow((ra1 - 1), 2) + pow((rb1 - 1), 2))); // v = 00, +1, +1//看卷积码寄存器状态
		gamma_pie_02 = exp(-1 / (2 * squared_sigma) * (pow((ra1 + 1), 2) + pow((rb1 + 1), 2))); // v = 11, -1, -1
		gamma_pie_10 = exp(-1 / (2 * squared_sigma) * (pow((ra1 + 1), 2) + pow((rb1 + 1), 2))); // v = 11, -1,-1

		gamma_pie_12 = exp(-1 / (2 * squared_sigma) * (pow((ra1 - 1), 2) + pow((rb1 - 1), 2))); // v = 00, +1, +1
		gamma_pie_21 = exp(-1 / (2 * squared_sigma) * (pow((ra1 + 1), 2) + pow((rb1 - 1), 2))); // v = 10, -1, +1

		gamma_pie_23 = exp(-1 / (2 * squared_sigma) * (pow((ra1 - 1), 2) + pow((rb1 + 1), 2))); // v = 01, +1,-1
		gamma_pie_31 = exp(-1 / (2 * squared_sigma) * (pow((ra1 - 1), 2) + pow((rb1 + 1), 2))); // v = 01, +1, -1
		gamma_pie_33 = exp(-1 / (2 * squared_sigma) * (pow((ra1 + 1), 2) + pow((rb1 - 1), 2))); // v = 10, -1, +1
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
}

void decoder()
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
	}
}