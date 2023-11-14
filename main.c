/***************************************************
Channel Coding Course Work: conolutional codes
This program template has given the message generator, BPSK modulation, AWGN channel model and BPSK demodulation,
you should first determine the encoder structure, then define the message and codeword length, generate the state table, write the convolutional encoder and decoder.
 
If you have any question, please contact me via e-mail: yanglj39@mail2.sysu.edu.cn
***************************************************/


/**
修改说明：
1. 原代码中 state_num 更改为 reg_num，为寄存器数量；state_num 为寄存器组的状态数量，state_num = 2^reg_num
*/

// #define  _CRT_SECURE_NO_WARNINGS // vs取消scanf报错
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "utils.h"

void statetable();
void encoder();
void modulation();
void demodulation();
void channel();
void decoder();

#define message_length 20 //the length of message
#define codeword_length 40 //the length of codeword
float code_rate; // 码率

// channel coefficient
#define pi 3.1415926
#define INF 0xFFFFFFF
double N0, sgm; // 信道噪声

int reg_num;//the number of the register of encoder structure
int state_num;//the number of the state of encoder structure
int input_len; // 每次输入的比特数，在758中为1
int input_states;
int* state_table;//state table, the size should be defined yourself

int message[message_length], codeword[codeword_length];//message and codeword，原信息和编码之后的信息
int re_codeword[codeword_length];//the received codeword
int de_message[message_length];//the decoding message

double tx_symbol[codeword_length][2];//the transmitted symbols
double rx_symbol[codeword_length][2];//the received symbols

void main()
{
	code_rate = (float)message_length / (float)codeword_length;
	PARAMETER* p = get_essential_params(7, 5, 8);
	reg_num = p->reg_num;//the number of the register of encoder structure
	state_num = pow(2, reg_num);//the number of the state of encoder structure
	input_len = 1; // 每次输入的比特数，在758中为1
	input_states = pow(2, input_len);


	int i;
	float SNR, start, finish;
	long int bit_error, seq, seq_num;
	double BER;
	double progress;

	//generate state table
	statetable();

	//random seed
	srand((int)time(0));

	// //input the SNR and frame number
	// printf("\nEnter start SNR: ");
	// scanf("%f", &start);
	// printf("\nEnter finish SNR: ");
	// scanf("%f", &finish);
	// printf("\nPlease input the number of message: ");
	// scanf("%d", &seq_num);
	start = 10, finish = 20; // 起始和结束的SNR，浮点数，单位为dB
	int SNR_step = 1; // SNR步长
	seq_num = 5; // 仿真次数

	for (SNR = start; SNR <= finish; SNR+=SNR_step)
	{
		//channel noise
		N0 = (1.0 / code_rate) / pow(10.0, (float)(SNR) / 10.0); // TODO：为什么是1/code_rate
		sgm = sqrt(N0 / 2);
		
		bit_error = 0;

		for (seq = 1; seq<=seq_num; seq++)
		{
			//generate binary message randomly
			/****************
			Pay attention that message is appended by 0 whose number is equal to the state of encoder structure.
			****************/
			// 随机生成信源序列
			for (i = 0; i<message_length - reg_num; i++)
			{
				message[i] = rand() % 2;
			}
			// 信源序列补零
			for (i = message_length - reg_num; i<message_length; i++)
			{
				message[i] = 0;
			}

			//convolutional encoder
			encoder();

			//BPSK modulation
			modulation();

			//AWGN channel
			// channel(); // TODO: 记得改回来

			//BPSK demodulation, it's needed in hard-decision Viterbi decoder
			//demodulation();

			//convolutional decoder
			decoder();

			//calculate the number of bit error
			for (i = 0; i<message_length; i++)
			{
				if (message[i] != de_message[i])
					bit_error++;
			}

			progress = (double)(seq * 100) / (double)seq_num;

			//calculate the BER
			BER = (double)bit_error / (double)(message_length*seq);

			//print the intermediate result
			printf("Progress=%2.1f, SNR=%2.1f, Bit Errors=%2.1d, BER=%E\r", progress, SNR, bit_error, BER);
		}

		//calculate the BER
		BER = (double)bit_error / (double)(message_length*seq_num);

		printf("Progress=%2.1f, SNR=%2.1f, Bit Errors=%2.1d, BER=%E\n", progress, SNR, bit_error, BER);
		if (DEBUG_MODE) {
			printf("[DEBUG]:\n");
			printf("信源序列:");
			for (int i=0; i<message_length; i++){
				printf("%d ", message[i]);
			}
			printf("\n");
			printf("编码后序列:");
			for (int i=0; i<codeword_length; i++){
				printf("%d ", codeword[i]);
			}
			printf("\n\n");
		}
	}
	// system("pause");
}


void statetable()
{
	// 把state_table填满
	// 仿照课件，state_table有state*2行5列
	// 列：IN, Current State, Next State, Out

	// (7, 5, 8)的state_table
	// 00,01... 用 0,1.. 来表示
	// TODO：之后补好版本的
	state_table = malloc(sizeof(int) * state_num * input_states * 5);
	for (int state=0; state<state_num; state++) {
		for (int input=0; input<input_states; input++) {
			int *inputArray = int2array(input, 1, BIN);
			int *stateArray = int2array(state, 2, BIN);

			int line_num = state*input_states+input;
			state_table[5*line_num] = input; // IN
			state_table[5*line_num+1] = state; // Current State

			int *nextStateArray = malloc(sizeof(int) * 2);
			nextStateArray[0] = inputArray[0];
			nextStateArray[1] = stateArray[0];
			state_table[5*line_num+2] = array2int(nextStateArray, 2); // Next State

			int *outArray = malloc(sizeof(int) * 2);
			outArray[0] = inputArray[0] ^ stateArray[0] ^ stateArray[1];
			outArray[1] = inputArray[0] ^ stateArray[1];
			state_table[5*line_num+3] = array2int(outArray, 2); // Out

			state_table[5*line_num+4] = line_num+1; // Line Number

			free(inputArray);
			free(stateArray);
			free(nextStateArray);
			free(outArray);
		}
	}

	// print state_table
	if (DEBUG_MODE) {
		printf("[DEBUG]:\n");
		printf("state_table:\n");
		for (int i=0; i<state_num*input_states; i++) {
			for (int j=0; j<5; j++) {
				printf("%d ", state_table[i*5+j]);
			}
			printf("\n");
		}
		printf("\n");
	}
}

void encoder()
{
	//convolution encoder, the input is message[] and the output is codeword[]

	// 先写一个最简单的(7, 5, 8)
	// state: s0, s1
	// input: u
	// output: c1 = u + s0 + s1, c2 = u + s1
	int s0 = 0, s1 = 0;
	for (int ii = 0; ii < message_length; ii++){
		codeword[ii*2] = (message[ii] + s0 + s1) % 2;
		codeword[ii*2+1] = (message[ii] + s1) % 2;
		// 移位
		s1 = s0;
		s0 = message[ii];
	}

}

void modulation()
{
	//BPSK modulation
	int i;

	//0 is mapped to (1,0) and 1 is mapped tp (-1,0)
	for (i = 0; i<codeword_length; i++)
	{
		tx_symbol[i][0] = -1 * (2 * codeword[i] - 1);
		tx_symbol[i][1]=0;
	}
}
void channel()
{
	//AWGN channel
	int i, j;
	double u, r, g;

	for (i = 0; i<codeword_length; i++)
	{
		for (j = 0; j<2; j++)
		{
			u=(float)rand()/(float)RAND_MAX;
			if(u==1.0)
				u=0.999999;
			r=sgm*sqrt(2.0*log(1.0/(1.0-u)));

			u=(float)rand()/(float)RAND_MAX;
			if(u==1.0)
				u=0.999999;
			g=(float)r*cos(2*pi*u);

			rx_symbol[i][j]=tx_symbol[i][j]+g;
		}
	}
}
void demodulation()
{
	int i;
	double d1, d2;
	for (i = 0; i<codeword_length; i++)
	{
		d1 = (rx_symbol[i][0] - 1)*(rx_symbol[i][0] - 1) + rx_symbol[i][1] * rx_symbol[i][1];
		d2 = (rx_symbol[i][0] + 1)*(rx_symbol[i][0] + 1) + rx_symbol[i][1] * rx_symbol[i][1];
		if (d1<d2)
			re_codeword[i] = 0;
		else
			re_codeword[i] = 1;
	}
}

int get_;

void decoder()
{ 
	// input: codeword
	// output: de_message

	// Viterbi decode
	// 硬判决维特比译码

} 