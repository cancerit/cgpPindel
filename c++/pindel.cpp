/*
 * This File is part of Pindel; a program to locate genomic variation.
 * https://trac.nbic.nl/pindel/
 * This is a CGP variant of 2.0.
 *
 *   Copyright (C) 2011 Kai Ye
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <iostream>
#include <fstream>
#include <string>
#include <string.h>
#include <vector>
#include <stdio.h>
#include <cstdlib>
#include <list>
#include <iomanip>
#include <cmath>
#include <math.h>
#include <time.h>
#include <bits/basic_string.h>
#include <map>
#include <stdint.h>

using namespace std;

const unsigned SpacerBeforeAfter = 10000000;

short Before, After;

unsigned int CountIndels = 0;
const int alphs = 4;
const char alphabet[alphs] = {'A','C','G','T'};
unsigned long long int TheMax = 0;
const short MAX_MISMATCHES = 4;
float ExtraDistanceRate = 0.1;

double Const_Log_T = 0.0;
double Const_S = 0.0;
double LOG14 = log10(0.25);
double Const_I = 0.0;
const unsigned int BoxSize = 10000;
const double Min_Filter_Ratio = 0.5;
unsigned int SPACERSIZE = 1;
unsigned int OriginalNumRead = 0;
const string NonACGT = "$";
short MIN_Len_Match = 4;
unsigned int NumberOfSIsInstances = 0;
unsigned int NumberOfDeletionsInstances = 0;
unsigned int NumberOfDIInstances = 0;
unsigned int NumberOfInvInstances = 0;
unsigned int NumberOfTDInstances = 0;
short ReportLength = 80;
vector <string> VectorTag;
char Match[256];
char Match2N[256];
char Convert2RC[256];
char Convert2RC4N[256];
char Cap2LowArray[256];
bool FirstChr = true;
const double InsertSizeExtra = 2;
unsigned int CONS_Chr_Size;
const string N_str = "N";
const char N_char = 'N';
const char BD_char = 'A';
unsigned int DSizeArray[15];

string BreakDancerMask;

const char Plus = '+';
const char Minus = '-';
const char FirstCharReadName = '@';

unsigned int NumReadScanned = 0;
unsigned int NumReadInChr = 0;
unsigned int InChrPlus = 0;
unsigned int InChrMinus = 0;
unsigned int GetPlus = 0;
unsigned int GetMinus = 0;

//short MAX_SNP_ERROR = 2;
const short ADDITIONAL_MISMATCH = 1;

//short TOTAL_SNP_ERROR_CHECKED = MAX_SNP_ERROR + ADDITIONAL_MISMATCH + 1;
//short TOTAL_SNP_ERROR_CHECKED_Minus = MAX_SNP_ERROR + ADDITIONAL_MISMATCH;

//short MAX_ALLOWED_MISMATCHES = TOTAL_SNP_ERROR_CHECKED_Minus + 5;

// #########################################################
short Min_Perfect_Match_Around_BP = 3;                   //#
const short MIN_IndelSize_NT = 3;                        //#
const short MIN_IndelSize_Inversion = 3;                 //#
const float Seq_Error_Rate = 0.05;                             //#
//float Seq_Error_Rate_1 = 0.05;                           //#
//float Seq_Error_Rate_2 = 0.02;                           //#
//float Seq_Error_Rate_3 = 0.00;                           //#
unsigned int BalanceCutoff = 50;                         //#
//short RangeMaxSensivity = 9;       // 3                  //#
//short RangeMediumSensivity = 9;    // 5                  //#
//short RangeLowSensivity = 7;                           //#
const bool Analyze_TD_INV_LI_Others = false;                    //#
const unsigned int NumRead2ReportCutOff = 3;             //#
short MaxRangeIndex = 9;// 5 or 6 or 7 or maximum 8//#
const float MaximumAllowedMismatchRate = 0.1;            //#
const short Min_Num_Matched_Bases = 30;                  //#
// #########################################################
//const float Double_Seq_Error_Rate_Per_Side = Seq_Error_Rate_Per_Side * 2;
unsigned int Distance = 300;
// short MinClose = 8;//short(log((double)Distance)/log(4.0) + 0.8) + 3 + MAX_SNP_ERROR;//atoi(argv[1]);
// short MinFar_I = MinClose + 1;//atoi(argv[2]);
//cout << "For short insertion: " << MinClose << "\t" << MinFar_I << endl;
 short MinFar_D = 8;//atoi(argv[3]);
const short MaxDI = 30;
const short FirstBase = 1;

//struct Fragment {
//   unsigned int Start;
//   unsigned int End;
//};



struct UniquePoint {
	short LengthStr;
	unsigned int AbsLoc;
	char Direction; // forward reverse
	char Strand; // sense antisense
	short Mismatches;
};

const char SENSE = '+';
const char ANTISENSE = '-';
const char FORWARD = '+';
const char BACKWARD = '-';

// UP_Close[0].AbsLoc

struct SPLIT_READ {
	string FragName;
	string Name;
	string UnmatchedSeq;
	char MatchedD;
	int MatchedRelPos;
	short MS;
	short InsertSize;
	string Tag;
	vector <UniquePoint> UP_Close; // partial alignment of the unmapped reads close to the mapped read
	vector <UniquePoint> UP_Far;
	vector <UniquePoint> UP_Far_backup;
	short ReadLength;
	short ReadLengthMinus;
	short MAX_SNP_ERROR;// = (short)(Temp_One_Read.UnmatchedSeq.size() * Seq_Error_Rate);

	short TOTAL_SNP_ERROR_CHECKED;// = MAX_SNP_ERROR + ADDITIONAL_MISMATCH + 1;
	short TOTAL_SNP_ERROR_CHECKED_Minus;// = MAX_SNP_ERROR + ADDITIONAL_MISMATCH;
	short MinClose;
	short BP;
	int Left;
	int Right;
	int BPLeft;
	int BPRight;
	int IndelSize;
	bool OK;
	double score;
	string InsertedStr;
	string NT_str;
	short NT_size;
	//string NT_2str;
	//short NT_2size;
	bool Used;
	short CloseEndLength;
};

struct LI_Pos {
	unsigned Plus_Pos;
	unsigned Minus_Pos;
	vector <unsigned> Plus_Reads; // put index here
	vector <unsigned> Minus_Reads;
	bool WhetherReport;
};

struct Rest_Pos {
	char Strand;
	unsigned Pos;
	vector <unsigned> Pos_Reads; // put index here
};

struct Indel4output {
   unsigned int BPLeft;
   unsigned int BPRight;
   int IndelSize;
   unsigned int Start;
   unsigned int End;
   unsigned int RealStart;
   unsigned int RealEnd;
   short NT_size;
   bool WhetherReport;
   string IndelStr;
   short Support;
};

void ReadInOneChr(ifstream & inf_Seq, string & TheInput, const string & ChrName);
short ReadInRead(ifstream & inf_Seq, const string & CurrentFragName, const string & CurrentFrag, vector <SPLIT_READ> & Reads);

void GetCloseEnd(const string & CurrentChr, SPLIT_READ & Temp_One_Read);
void GetFarEnd_BothStrands(const string & CurrentChr, SPLIT_READ & Temp_One_Read, const short & RangeIndex);
void GetFarEnd_OtherStrand(const string & CurrentChr, SPLIT_READ & Temp_One_Read, const short & RangeIndex);
void GetFarEnd_SingleStrandDownStream(const string & CurrentChr, SPLIT_READ & Temp_One_Read, const short & RangeIndex);
void GetFarEnd_SingleStrandUpStream(const string & CurrentChr, SPLIT_READ & Temp_One_Read, const short & RangeIndex);
void GetFarEnd_SingleStrandDownStreamInsertions(const string & CurrentChr, SPLIT_READ & Temp_One_Read, const short & RangeIndex);
//void GetFarEnd(const string & CurrentChr, SPLIT_READ & Temp_One_Read, const int & Range, const string & BreakDancerMask);
void GetFarEnd(const string & CurrentChr, SPLIT_READ & Temp_One_Read, const int & start, const int & end);

string GetConsensusInsertedStr(const vector <SPLIT_READ> & Reads, const int & StartIndex, const int & EndIndex);
struct BreakDancer {
	string ChrName_A;
	string ChrName_B;
	string Type;
	int Size;
	int Score;
	unsigned S1;
	unsigned S2;
	unsigned S3;
	unsigned S4;
};

struct Region {
	unsigned start;
	unsigned end;
};

void OutputDeletions(const vector <SPLIT_READ> & Deletions,
                     const string & TheInput,
                     const unsigned int & C_S,
                     const unsigned int & C_E,
		               const unsigned int & RealStart,
		               const unsigned int & RealEnd,
                     ofstream & DeletionOutf);

void OutputTDs(const vector <SPLIT_READ> & TDs,
                     const string & TheInput,
                     const unsigned int & C_S,
                     const unsigned int & C_E,
		               const unsigned int & RealStart,
		               const unsigned int & RealEnd,
                     ofstream & TDOutf);


void OutputInversions(const vector <SPLIT_READ> & Inv,
                     const string & TheInput,
                     const unsigned int & C_S,
                     const unsigned int & C_E,
		               const unsigned int & RealStart,
		               const unsigned int & RealEnd,
                     ofstream & InvOutf);

void OutputSIs(const vector <SPLIT_READ> & SIs,
               const string & TheInput,
               const unsigned int & C_S,
               const unsigned int & C_E,
	            const unsigned int & RealStart,
	            const unsigned int & RealEnd,
               ofstream & SIsOutf);



//void SortAndPutLI(vector <SPLIT_READ> & Input, vector <int> & Output, const string & CurrentFragName, const string & chr, const char & Direction, ofstream & LI_SR_Outf);
//void ProcessLIs(vector <SPLIT_READ> & LIs, ofstream & LIoutputfile);

short CompareTwoReads(const SPLIT_READ & First, const SPLIT_READ & Second);
vector <string> ReverseComplement(const vector <string> & input);
string Reverse(const string & InputPattern);
string ReverseComplement(const string & InputPattern);
string Cap2Low(const string & input);
void CheckLeft_Close(const SPLIT_READ & OneRead,
							const string & TheInput,
               const string & CurrentReadSeq,
               const vector <unsigned int> Left_PD[],
               const short & BP_Left_Start,
               const short & BP_Left_End,
               const short & CurrentLength,
               vector <UniquePoint> & LeftUP);

void CheckRight_Close(const SPLIT_READ & OneRead,
							 const string & TheInput,
                const string & CurrentReadSeq,
                const vector <unsigned int> Right_PD[],
                const short & BP_Right_Start,
                const short & BP_Right_End,
                const short & CurrentPos,
                vector <UniquePoint> & RightUP);

void CheckLeft_Far(const SPLIT_READ & OneRead,
						 const string & TheInput,
							const string & CurrentReadSeq,
							const vector <unsigned int> Left_PD[],
							const short & BP_Left_Start,
							const short & BP_Left_End,
							const short & CurrentLength,
							vector <UniquePoint> & LeftUP);

void CheckRight_Far(const SPLIT_READ & OneRead,
						  const string & TheInput,
					  const string & CurrentReadSeq,
					  const vector <unsigned int> Right_PD[],
					  const short & BP_Right_Start,
					  const short & BP_Right_End,
					  const short & CurrentPos,
					  vector <UniquePoint> & RightUP);

void CheckBoth(const SPLIT_READ & OneRead,
					const string & TheInput,
               const string & CurrentReadSeq,
               const vector <unsigned int> PD_Plus[],
					const vector <unsigned int> PD_Minus[],
               const short & BP_Start,
               const short & BP_End,
               const short & CurrentLength,
               vector <UniquePoint> & UP);
bool CheckMismatches(const string & TheInput,
                     const string & CurrentReadSeq,
                     //const unsigned int & Start,
							const UniquePoint & UP);

//void ProcessLIs(vector <SPLIT_READ> & LIs, ofstream & LIoutputfile);
void SortOutputSI(const unsigned & NumBoxes, const string & CurrentChr, vector <SPLIT_READ> & AllReads, vector <unsigned> SIs[], ofstream & SIsOutf);
void SortOutputD(const unsigned & NumBoxes, const string & CurrentChr, vector <SPLIT_READ> & AllReads, vector <unsigned> Deletions[], ofstream & DeletionOutf);
void SortOutputTD(const unsigned & NumBoxes, const string & CurrentChr, vector <SPLIT_READ> & AllReads, vector <unsigned> TDs[], ofstream & TDOutf);
void SortOutputTD_NT(const unsigned & NumBoxes, const string & CurrentChr, vector <SPLIT_READ> & AllReads, vector <unsigned> TDs[], ofstream & TDOutf);
void SortOutputInv(const unsigned & NumBoxes, const string & CurrentChr, vector <SPLIT_READ> & AllReads, vector <unsigned> Inv[], ofstream & InvOutf);
void SortOutputInv_NT(const unsigned & NumBoxes, const string & CurrentChr, vector <SPLIT_READ> & AllReads, vector <unsigned> Inv[], ofstream & InvOutf);
void SortOutputDI(const unsigned & NumBoxes, const string & CurrentChr, vector <SPLIT_READ> & AllReads, vector <unsigned> DI[], ofstream & DIOutf);
void SortOutputLI(const string & CurrentChr, vector <SPLIT_READ> & AllReads, ofstream & Outf_LI);
void SortOutputRest(const string & CurrentChr, vector <SPLIT_READ> & AllReads, ofstream & Outf_Rest);
bool NotInVector(const string & OneTag, const vector <string> & VectorTag);
void GetIndelTypeAndRealStart(const string & TheInput, const unsigned int & BPLeft,
                              const unsigned int & IndelSize, const string & IndelStr,
                              string & IndelType, unsigned int & RealStart, const bool & WhetherD);
void GetRealStart4Deletion(const string & TheInput, unsigned int & RealStart, unsigned int & RealEnd);
void GetRealStart4Insertion(const string & TheInput, string & InsertedStr, unsigned int & RealStart, unsigned int & RealEnd);
bool ReportEvent(const vector <SPLIT_READ> & Deletions, const unsigned int & Pre_S, const unsigned int & Pre_E);

void CleanUniquePoints (vector <UniquePoint> & Input_UP);

short WhetherRemoveDuplicates;

string TempLine_DB_OK;

vector <Region> Merge(const vector <Region> & AllRegions);

short CompareTwoString(const string & Str_A, const string & Str_B);

int main (int argc, char *argv[]) {

	if (NumRead2ReportCutOff == 1) BalanceCutoff = 3000000000;

	// ################### module 1: define input and output #####################
  if (argc != 7) {
    cout << "\nWelcome to Pindel, developed by Kai Ye, k.ye@lumc.nl\n\n"
         << "5 parameters are required here:\n"
         << "1. Input: the reference genome sequences in fasta format;\n"
         << "2. Input: the unmapped reads in a modified fastq format;\n"
	      << "   Better to use bam2pindel.pl to convert BAM files to Pindel input.\n"
	      << "   If the perl script fails for some reasons, please use the provided\n"
	      << "   sam2pindel.cpp to extract reads from sam files.\n"
	      << "   Compile cpp file first: g++ sam2pindel.cpp -o sam2pindel -O3\n"
	      << "3. Output folder\n"
			<< "4. Which chr/fragment\n"
	      << "   Pindel will process reads for one chr each time\n"
	      << "   ChrName must be the same as in reference sequence and in read file\n"
			<< "5. BreakDancer result:\n"
	      << "   ChrA LocA stringA ChrB LocB stringB others\n"
	      << "   If you don't have BreakDancer result, please provide an empty file here.\n"
	      << "6. Maximum event size index. 5 is recommended\n" // MaxRangeIndex = 5
	      << "   2:          128\n"
	      << "   3:          512\n"
	      << "   4:        2,048\n"
	      << "   5:        8,092\n"
    	   << "   6:       32,368\n"
	      << "   7:      129,472\n"
	      << "   8:      517,888\n"
	      << "   9:    2,071,552\n"
	      << "   10:   8,286,208\n"
	      << "   11:  33,144,832\n"
	      << "   12: 132,579,328\n"
	      //<< "   13: 530,317,312\n"

         << endl;
    return 1;
  }

  // #################################################################
	MaxRangeIndex = atoi(argv[6]);
	if (MaxRangeIndex <= 2) {
		cout << "Please set Maximum event size larger than 2" << endl;
		return 0;
	}
	if (MaxRangeIndex > 12) {
		cout << "Please set Maximum event size <= 12" << endl;
		return 0;
	}
   if (MaxRangeIndex > 5 && MaxRangeIndex <= 8) cout << "Pindel may be slow if you have too much data." << endl;
	if (MaxRangeIndex > 8) cout << "Pindel may be VERY slow if you have too much data." << endl;

	const string WhichChr = argv[4];

	ifstream  inf_Seq(argv[1]);   // input file name

	 string inputFilename = argv[2];
	inputFilename += "_" + WhichChr;
	bool WithFolder = false;
	int StartOfFileName = 0;
	for (int i = inputFilename.size(); i >= 0; i--) {
		if (inputFilename[i] == '/') {
			StartOfFileName = i;
			WithFolder = true;
			break;
		}
	}

	if (WithFolder) {
		inputFilename = inputFilename.substr(StartOfFileName + 1, inputFilename.size() - 1 - StartOfFileName);
	}
	//cout << inputFilename << endl;


	 string OutputFolder = argv[3];

	 OutputFolder += "/";

    string SIOutputFilename = OutputFolder + inputFilename + "_SI";      // output file name
    //strcpy(SIOutputFilename, (OutputFolder + inputFilename + "_SI").c_str());
    ofstream SIoutputfile_test(SIOutputFilename.c_str());
	 if (!SIoutputfile_test) {
		 cout << "Sorry, cannot write to the file: " << SIOutputFilename << endl;
		 return 1;
	 }
	 SIoutputfile_test.close();

    //char DeletionOutputFilename[10000];      // output file name
	 string DeletionOutputFilename = OutputFolder + inputFilename + "_D";
    //strcpy(DeletionOutputFilename, (OutputFolder + inputFilename + "_D").c_str());
    ofstream DeletionOutf_test(DeletionOutputFilename.c_str());
	 if (!DeletionOutf_test) {
		 cout << "Sorry, cannot write to the file: " << DeletionOutputFilename << endl;
		 return 1;
	 }
	 DeletionOutf_test.close();

	string TDOutputFilename = OutputFolder + inputFilename + "_TD";
	//strcpy(DeletionInsertinOutputFilename, (OutputFolder + inputFilename + "_DI").c_str());
	ofstream TDOutf_test(TDOutputFilename.c_str());
	if (!TDOutf_test) {
		cout << "Sorry, cannot write to the file: " << TDOutputFilename << endl;
		return 1;
	}
	TDOutf_test.close();

	//char InversionOutputFilename[10000];      // output file name
	string InversionOutputFilename = OutputFolder + inputFilename + "_INV";
	//strcpy(InversionOutputFilename, (OutputFolder + inputFilename + "_INV").c_str());
	ofstream InversionOutf_test(InversionOutputFilename.c_str());
	if (!InversionOutf_test) {
		cout << "Sorry, cannot write to the file: " << InversionOutputFilename << endl;
		return 1;
	}
	InversionOutf_test.close();

	//char LargeInsertionOutputFilename[10000];      // output file name
	string LargeInsertionOutputFilename = OutputFolder + inputFilename + "_LI";
	//strcpy(LargeInsertionOutputFilename, (OutputFolder + inputFilename + "_LI").c_str());
	ofstream LargeInsertionOutf_test(LargeInsertionOutputFilename.c_str());
	if (!LargeInsertionOutf_test) {
		cout << "Sorry, cannot write to the file: " << LargeInsertionOutputFilename << endl;
		return 1;
	}
	LargeInsertionOutf_test.close();

	//char RestOutputFilename[10000];      // output file name
	string RestOutputFilename = OutputFolder + inputFilename + "_BP";
	//strcpy(RestOutputFilename, (OutputFolder + inputFilename + "_BP").c_str());
	ofstream RestOutf_test(RestOutputFilename.c_str());
	if (!RestOutf_test) {
		cout << "Sorry, cannot write to the file: " << RestOutputFilename << endl;
		return 1;
	}
	RestOutf_test.close();



	//WhetherRemoveDuplicates = atoi(argv[6]);
	int Count_SI = 0;
	int Count_D = 0;
	int Count_DI = 0;
	int Count_TD = 0;
	int Count_TD_NT = 0;
	int Count_Inv = 0;
	int Count_Inv_NT = 0;
	int Count_D_Plus = 0;
	int Count_D_Minus = 0;
	int Count_DI_Plus = 0;
	int Count_DI_Minus = 0;
	int Count_TD_Plus = 0;
	int Count_TD_Minus = 0;
	int Count_TD_NT_Plus = 0;
	int Count_TD_NT_Minus = 0;
	int Count_Inv_Plus = 0;
	int Count_Inv_Minus = 0;
	int Count_Inv_NT_Plus = 0;
	int Count_Inv_NT_Minus = 0;
	int Count_SI_Plus = 0;
	int Count_SI_Minus = 0;
	//int Entering_D_Plus = 0;
	//int Entering_D_Minus = 0;
	//int Plus_Sum_Left = 0;
	//int Plus_Sum_Right = 0;
	//int Minus_Sum_Left = 0;
	//int Minus_Sum_Right = 0;

	Match[(short)'A'] = 'A';
	Match[(short)'C'] = 'C';
	Match[(short)'G'] = 'G';
	Match[(short)'T'] = 'T';
	Match[(short)'N'] = 'X';
	Match[(short)'$'] = '$';
	Match2N[(short)'A'] = 'N';
	Match2N[(short)'C'] = 'N';
	Match2N[(short)'G'] = 'N';
	Match2N[(short)'T'] = 'N';
	Match2N[(short)'N'] = 'X';
	Match2N[(short)'$'] = '$';
	Convert2RC[(short)'A'] = 'T';
	Convert2RC[(short)'C'] = 'G';
	Convert2RC[(short)'G'] = 'C';
	Convert2RC[(short)'T'] = 'A';
	Convert2RC[(short)'N'] = 'X';
	Convert2RC[(short)'$'] = '$';
	Convert2RC4N[(short)'A'] = 'T';
	Convert2RC4N[(short)'C'] = 'G';
	Convert2RC4N[(short)'G'] = 'C';
	Convert2RC4N[(short)'T'] = 'A';
	Convert2RC4N[(short)'N'] = 'N';
	Cap2LowArray[(short)'A'] = 'a';
	Cap2LowArray[(short)'C'] = 'c';
	Cap2LowArray[(short)'G'] = 'g';
	Cap2LowArray[(short)'T'] = 't';
	Cap2LowArray[(short)'N'] = 'n';
	Cap2LowArray[(short)'$'] = 'n';


	time_t Time_Load_S, Time_Load_E, Time_Mine_E, Time_Sort_E;//, Time_End;
	Time_Load_S = time(NULL);
	unsigned int AllLoadings = 0;
	unsigned int AllMinings = 0;
	unsigned int AllSortReport = 0;


    string Spacer = "";
    for (int i = 0; i < SpacerBeforeAfter; i++)
       Spacer += "N";
    //cout << Distance << endl;
	 Distance = 300;


	//DSizeArray[0]  = 0;
	DSizeArray[0]  =         25;
	DSizeArray[1]  = 128;
	DSizeArray[2]  = DSizeArray[1] * 4; // 512
	DSizeArray[3]  = DSizeArray[2] * 4; // 2048
	DSizeArray[4]  = DSizeArray[3] * 4; // 8092
	DSizeArray[5]  = DSizeArray[4] * 4; // 32368
	DSizeArray[6]  = DSizeArray[5] * 4; // 129472
	DSizeArray[7]  = DSizeArray[6] * 4; // 517888
	DSizeArray[8]  = DSizeArray[7] * 4;

	DSizeArray[9]  = DSizeArray[8] * 4;
	DSizeArray[10] = DSizeArray[9] * 4;
	DSizeArray[11] = DSizeArray[10] * 4;
	DSizeArray[12] = DSizeArray[11] * 4;
	DSizeArray[13] = DSizeArray[12] * 4;
	DSizeArray[14] = DSizeArray[13] * 4;

    //unsigned int DSizeExtra = 100;
    //unsigned int D_SIZE = 100;

	string TempLie_BD;

    string CurrentChr;
    char TempStartChrChar;
    //short FragID;

	char FirstSharpChar;

    unsigned int EndOfFragment;// = CurrentChr.size() - Spacer;
    unsigned int StartOfFragment;// = CurrentChr.size() - Spacer;

	//short BP_Left_Start;
	//short BP_Left_End;
	//short BP_Right_Start;
	//short BP_Right_End;
    //unsigned int LeftStart, LeftEnd, RightStart, RightEnd;
    //bool ShortInsertionOK;
    //bool DeletionOK;
    //unsigned int DISTANCE;
    //short ReadLengthMinus;
    //short ReadLength;

    //char LeftChar, RightChar;
    unsigned int Num_Left;
    //short CurrentChrIndex = 0;


    // ################### module 2: read input #####################

	string TempLine_BD;

    //while (inf_Seq >> CurrentFragName)

       ReadInOneChr(inf_Seq, CurrentChr, WhichChr);
	    if (CurrentChr.empty()) {
			 cout << "Cannot find the requested chr." << endl;
			 return 1;
		 }
	    CONS_Chr_Size = CurrentChr.size() - 2 * SpacerBeforeAfter;
       unsigned NumBoxes =(unsigned) (CurrentChr.size() / BoxSize) + 1;
		 cout << NumBoxes << "\t" << BoxSize << endl;

      // cout << "NumBoxes: " << NumBoxes << endl;
       vector <unsigned> SIs[NumBoxes];
       //vector <SPLIT_READ> LIs[NumBoxes];
       vector <unsigned> Deletions[NumBoxes];
		 vector <unsigned> TD[NumBoxes];
		 vector <unsigned> TD_NT[NumBoxes];
       vector <unsigned> DI[NumBoxes];
		 vector <unsigned> Inv[NumBoxes];
		 vector <unsigned> Inv_NT[NumBoxes];

       EndOfFragment = CurrentChr.size() - SpacerBeforeAfter;
       StartOfFragment = SpacerBeforeAfter;
       //cout << StartOfFragment << "\t" << EndOfFragment << endl;
       vector <SPLIT_READ> InputReads, Reads;
       ifstream  inf_ReadsSeq(argv[2]);   // input file name

		 ifstream  inf_BP_test(argv[5]);   // input file name
		 ifstream  inf_BP(argv[5]);   // input file name
		 vector <BreakDancer> All_BD_events;
		 BreakDancer Temp_BD_event;
		 All_BD_events.push_back(Temp_BD_event);

	    //short * BD_INDEX;// = new short[CurrentChr.size()];


		 //if (WhichChr == CurrentFragName)
		 //{
			 while (inf_BP_test >> FirstSharpChar) {
				 if (FirstSharpChar == '#') {
					 getline(inf_BP_test, TempLine_BD);
					 getline(inf_BP, TempLine_BD);
				 }
				 else {
					 getline(inf_BP_test, TempLine_BD);
					 inf_BP >> Temp_BD_event.ChrName_A >> Temp_BD_event.S1 >> TempLine_BD
					        >> Temp_BD_event.ChrName_B >> Temp_BD_event.S3 >> TempLine_BD;
					 //>> Temp_BD_event.Type >> Temp_BD_event.Size;
					 getline(inf_BP, TempLine_BD);

					 //<< Temp_BD_event.Size << endl;
					 //if (Temp_BD_event.Type == "DEL" )
					 if (WhichChr == Temp_BD_event.ChrName_A)
					 {
                   if (Temp_BD_event.S1 + 200 > CONS_Chr_Size) Temp_BD_event.S2 = CONS_Chr_Size - 1;
						 else Temp_BD_event.S2 = Temp_BD_event.S1 + 200;
						 if (Temp_BD_event.S1 > 200) Temp_BD_event.S1 = Temp_BD_event.S1 - 200;
						 else Temp_BD_event.S1 = 1;
					 }
					 if (WhichChr == Temp_BD_event.ChrName_B)
					 {
                   if (Temp_BD_event.S3 + 200 > CONS_Chr_Size) Temp_BD_event.S4 = CONS_Chr_Size - 1;
						 else Temp_BD_event.S4 = Temp_BD_event.S3 + 200;
						 if (Temp_BD_event.S3 > 200) Temp_BD_event.S3 = Temp_BD_event.S3 - 200;
						 else Temp_BD_event.S3 = 1;
					 }
					 Temp_BD_event.S1 += SpacerBeforeAfter;
					 Temp_BD_event.S2 += SpacerBeforeAfter;
					 Temp_BD_event.S3 += SpacerBeforeAfter;
					 Temp_BD_event.S4 += SpacerBeforeAfter;
					 if (WhichChr == Temp_BD_event.ChrName_A && WhichChr == Temp_BD_event.ChrName_B) {
						 //cout << WhichChr << " " << Temp_BD_event.S1 << " " << Temp_BD_event.S2 << " " << Temp_BD_event.S3 << " " << Temp_BD_event.S4 << endl;
						 All_BD_events.push_back(Temp_BD_event);
					 }
				 }
			 }
			 cout << "BreakDancer events: " << All_BD_events.size() - 1 << endl;

		 //}

		 short ReturnFromReadingReads;
		 int GetPlus = 0;
		 int GetMinus = 0;
		 //if (WhichChr == CurrentFragName)
		 {
          ReturnFromReadingReads = 0;
          ReturnFromReadingReads = ReadInRead(inf_ReadsSeq, WhichChr, CurrentChr, Reads);
          if (ReturnFromReadingReads == 1) {
             cout << "malformed record detected!" << endl;
	          return 1;
          }
			 else if (Reads.size() == 0) return 0;
       }

		 Time_Mine_E = time(NULL);

		 if (Reads.size())
		    cout << "There are " << Reads.size() << " reads for this chromosome." <<endl;
       else {
          cout << "There are no reads for this chromosome." <<endl;
          return 0;
       }
       Num_Left = Reads.size();
       Const_Log_T = log10((double)Num_Left);
		 Time_Load_E = time(NULL);
		 int CountFarEnd, CountFarEndPlus, CountFarEndMinus;

		 // ################### module 3: search breakpoints #####################
		 if (All_BD_events.size() > 1) {
			 cout << "Searching additional breakpoints by adding BreakDancer results" << endl;
			    short * BD_INDEX = new short[CurrentChr.size()];
				 for (unsigned i = 0; i < CurrentChr.size(); i++) BD_INDEX[i] = 0;
				 for (unsigned i = 1; i < All_BD_events.size(); i++) {
					 //cout << i << endl;
					 for (unsigned j = All_BD_events[i].S1; j < All_BD_events[i].S2; j++) BD_INDEX[j] = i;
					 for (unsigned j = All_BD_events[i].S3; j < All_BD_events[i].S4; j++) BD_INDEX[j] = i * (-1);
				 }
				 int BD_Plus = 0;
				 int BD_Minus = 0;
				 for (unsigned i = 0; i < CurrentChr.size(); i++) {
					 if (BD_INDEX[i] > 0) BD_Plus++;
					 else if (BD_INDEX[i] < 0) BD_Minus++;
				 }
				 cout << BD_Plus << "\t" << BD_Minus << endl;

			 //ADDITIONAL_MISMATCH = 2;
			 //Seq_Error_Rate = 0.05;
			 //int BD_event_Index = 0;
			 CountFarEnd = 0;
			 CountFarEndMinus = 0;
			 CountFarEndPlus = 0;
			 int Start_pos, End_pos;
			 for (unsigned ReadIndex = 0; ReadIndex < Reads.size(); ReadIndex++) {
				 if (!Reads[ReadIndex].UP_Far.empty()) {
					 //CountFarEnd++;
					 continue;
				 }
				 int BD_event_Index = BD_INDEX[Reads[ReadIndex].UP_Close[0].AbsLoc];
				 if (BD_event_Index == 0) continue;
				 //All_BD_events[i].S4; j++) BD_INDEX[j]
				 //MAX_SNP_ERROR = (short)(Reads[ReadIndex].UnmatchedSeq.size() * Seq_Error_Rate + 1.0);

				 //TOTAL_SNP_ERROR_CHECKED = MAX_SNP_ERROR + ADDITIONAL_MISMATCH + 1;
				 //TOTAL_SNP_ERROR_CHECKED_Minus = MAX_SNP_ERROR + ADDITIONAL_MISMATCH;
				 //MinClose = short(log((double)(2 * Reads[ReadIndex].InsertSize)/log(4.0) + 0.8)) + 3;//atoi(argv[1]);
				 //MinFar_I = MinClose + 1;//atoi(argv[2]);
				 //cout << "For short insertion: " << MinClose << "\t" << MinFar_I << endl;
				 //MinFar_D = 8;//atoi(argv[3]);
				 if (BD_event_Index > 0) {
					 //cout << ReadIndex << "\t" << BD_event_Index << "\t" << All_BD_events[BD_event_Index].S1 << "\t" << All_BD_events[BD_event_Index].S2 << "\t" << All_BD_events[BD_event_Index].S3 << "\t" << All_BD_events[BD_event_Index].S4 << endl;
					 Start_pos = All_BD_events[BD_event_Index].S3 - Reads[ReadIndex].ReadLength;
					 End_pos = All_BD_events[BD_event_Index].S4 + Reads[ReadIndex].ReadLength;
					 //cout << Start_pos << "\t" << End_pos << endl;
					 GetFarEnd(CurrentChr, Reads[ReadIndex], Start_pos, End_pos);
				 }
				 else { // < 0
					 //cout << ReadIndex << "\t" << BD_event_Index << "\t" << All_BD_events[BD_event_Index * (-1)].S1 << "\t" << All_BD_events[BD_event_Index * (-1)].S2 << "\t" << All_BD_events[BD_event_Index * (-1)].S3 << "\t" << All_BD_events[BD_event_Index * (-1)].S4 << endl;
					 Start_pos = All_BD_events[BD_event_Index * (-1)].S1 - Reads[ReadIndex].ReadLength;
					 End_pos = All_BD_events[BD_event_Index * (-1)].S2 + Reads[ReadIndex].ReadLength;
					 //cout << Start_pos << "\t" << End_pos << endl;
					 GetFarEnd(CurrentChr, Reads[ReadIndex], Start_pos, End_pos);
				 }

				 if (!Reads[ReadIndex].UP_Far.empty()) {
					 //cout << ReadIndex << "\t" << Reads[ReadIndex].UP_Close[0].AbsLoc << "\t" << Reads[ReadIndex].UP_Far[0].AbsLoc << endl;
					 if (Reads[ReadIndex].UP_Far[Reads[ReadIndex].UP_Far.size() - 1].LengthStr + Reads[ReadIndex].CloseEndLength < Reads[ReadIndex].ReadLength) {
						 if (Reads[ReadIndex].UP_Far_backup.size()) {
							 if (Reads[ReadIndex].UP_Far_backup[Reads[ReadIndex].UP_Far_backup.size() - 1].LengthStr <Reads[ReadIndex].UP_Far[Reads[ReadIndex].UP_Far.size() - 1].LengthStr) {
								 Reads[ReadIndex].UP_Far_backup = Reads[ReadIndex].UP_Far;
								 Reads[ReadIndex].UP_Far.clear();
							 }
							 else Reads[ReadIndex].UP_Far.clear();
						 }
						 else {
							 Reads[ReadIndex].UP_Far_backup = Reads[ReadIndex].UP_Far;
							 Reads[ReadIndex].UP_Far.clear();
						 }
					 }
					 else {
						 CountFarEnd++;
						 if (Reads[ReadIndex].MatchedD == Plus) CountFarEndPlus++;
						 else CountFarEndMinus++;
					 }
				 }
			 }
			 cout << "\tNumber of reads with far end mapped: " << CountFarEnd << "\t"
			 //<< "\tTotal number of reads:" << Reads.size() << "\n"
			 //<< CountFarEnd * 100.0 / Reads.size() << " %\n"
			 << "Far+: " << CountFarEndPlus << "\tFar-: " << CountFarEndMinus << endl;
		    delete[] BD_INDEX;
		 }

		 cout << "Searching breakpoints of deletion events" << endl;
		 for (short RangeIndex = 1; RangeIndex < MaxRangeIndex; RangeIndex++) {

			 CountFarEnd = 0;
			 CountFarEndMinus = 0;
			 CountFarEndPlus = 0;
			 for (unsigned ReadIndex = 0; ReadIndex < Reads.size(); ReadIndex++) {
				 if (!Reads[ReadIndex].UP_Far.empty()) {
					 //CountFarEnd++;
					 continue;
				 }

				 //MinClose = short(log((double)(2 * Reads[ReadIndex].InsertSize + DSizeArray[RangeIndex]))/log(4.0) + 0.8) + 3;//atoi(argv[1]);
				 //MinFar_I = MinClose + 1;//atoi(argv[2]);
				 //cout << "For short insertion: " << MinClose << "\t" << MinFar_I << endl;
				 //MinFar_D = 8;//atoi(argv[3]);
				 //if (RangeIndex <= 5)
				    GetFarEnd_SingleStrandDownStream(CurrentChr, Reads[ReadIndex], RangeIndex);
				 if (!Reads[ReadIndex].UP_Far.empty()) {
					 if (Reads[ReadIndex].UP_Far[Reads[ReadIndex].UP_Far.size() - 1].LengthStr + Reads[ReadIndex].CloseEndLength < Reads[ReadIndex].ReadLength) {
						 if (Reads[ReadIndex].UP_Far_backup.size()) {
							 if (Reads[ReadIndex].UP_Far_backup[Reads[ReadIndex].UP_Far_backup.size() - 1].LengthStr <Reads[ReadIndex].UP_Far[Reads[ReadIndex].UP_Far.size() - 1].LengthStr) {
								 Reads[ReadIndex].UP_Far_backup = Reads[ReadIndex].UP_Far;
								 Reads[ReadIndex].UP_Far.clear();
							 }
							 else Reads[ReadIndex].UP_Far.clear();
						 }
						 else {
							 Reads[ReadIndex].UP_Far_backup = Reads[ReadIndex].UP_Far;
							 Reads[ReadIndex].UP_Far.clear();
						 }
					 }
					 else {
						 CountFarEnd++;
						 if (Reads[ReadIndex].MatchedD == Plus) CountFarEndPlus++;
						 else CountFarEndMinus++;
					 }
				 }
			 }

			 cout << RangeIndex << "\tNumber of reads with far end mapped: " << CountFarEnd << "\t"
			      //<< "\tTotal number of reads:" << Reads.size() << "\n"
			      //<< CountFarEnd * 100.0 / Reads.size() << " %\n"
			      << "Far+: " << CountFarEndPlus << "\tFar-: " << CountFarEndMinus << endl;
		 }

		 cout << "Searching breakpoints of SI events" << endl;
		 for (short RangeIndex = 1; RangeIndex < 2; RangeIndex++) {
			 CountFarEnd = 0;
			 CountFarEndMinus = 0;
			 CountFarEndPlus = 0;

			 for (unsigned ReadIndex = 0; ReadIndex < Reads.size(); ReadIndex++) {
				 if (!Reads[ReadIndex].UP_Far.empty()) {
					 //CountFarEnd++;
					 continue;
				 }
				 {
					 //MAX_SNP_ERROR = (short)(Reads[ReadIndex].UnmatchedSeq.size() * Seq_Error_Rate + 1.0);

					 //TOTAL_SNP_ERROR_CHECKED = MAX_SNP_ERROR + ADDITIONAL_MISMATCH + 1;
					 //TOTAL_SNP_ERROR_CHECKED_Minus = MAX_SNP_ERROR + ADDITIONAL_MISMATCH;
				 }
				 //MinClose = short(log((double)(2 * Reads[ReadIndex].InsertSize + DSizeArray[RangeIndex]))/log(4.0) + 0.8) + 3;//atoi(argv[1]);
				 //MinFar_I = MinClose + 1;//atoi(argv[2]);
				 //cout << "For short insertion: " << MinClose << "\t" << MinFar_I << endl;
				 //MinFar_D = 8;//atoi(argv[3]);
				 //if (RangeIndex <= 5)
				 GetFarEnd_SingleStrandDownStreamInsertions(CurrentChr, Reads[ReadIndex], RangeIndex);
				 if (!Reads[ReadIndex].UP_Far.empty()) {
				    CountFarEnd++;
					 if (Reads[ReadIndex].MatchedD == Plus) CountFarEndPlus++;
					 else CountFarEndMinus++;
				 }
			 }

			 cout << RangeIndex << "\tNumber of reads with far end mapped: " << CountFarEnd << "\t"
			 //<< "\tTotal number of reads:" << Reads.size() << "\n"
			 //<< CountFarEnd * 100.0 / Reads.size() << " %\n"
			 << "Far+: " << CountFarEndPlus << "\tFar-: " << CountFarEndMinus << endl;
		 }

if (Analyze_TD_INV_LI_Others) {


		 cout << "Searching breakpoints of tandem duplication events" << endl;
		 //int CountFarEnd, CountFarEndPlus, CountFarEndMinus;
		 for (short RangeIndex = 1; RangeIndex < MaxRangeIndex; RangeIndex++) {

			 CountFarEnd = 0;
			 CountFarEndMinus = 0;
			 CountFarEndPlus = 0;
			 for (unsigned ReadIndex = 0; ReadIndex < Reads.size(); ReadIndex++) {
				 if (!Reads[ReadIndex].UP_Far.empty()) {
					 //CountFarEnd++;
					 continue;
				 }

				 //MinClose = short(log((double)(2 * Reads[ReadIndex].InsertSize + DSizeArray[RangeIndex]))/log(4.0) + 0.8) + 3;//atoi(argv[1]);
				 //MinFar_I = MinClose + 1;//atoi(argv[2]);
				 //cout << "For short insertion: " << MinClose << "\t" << MinFar_I << endl;
				 //MinFar_D = 8;//atoi(argv[3]);
				 //if (RangeIndex <= 5)
				 GetFarEnd_SingleStrandUpStream(CurrentChr, Reads[ReadIndex], RangeIndex);
				 if (!Reads[ReadIndex].UP_Far.empty()) {
					 if (Reads[ReadIndex].UP_Far[Reads[ReadIndex].UP_Far.size() - 1].LengthStr + Reads[ReadIndex].CloseEndLength >= Reads[ReadIndex].ReadLength) {
						 if (Reads[ReadIndex].MatchedD == Plus) {
							 if (Reads[ReadIndex].UP_Close[0].AbsLoc < Reads[ReadIndex].ReadLength + Reads[ReadIndex].UP_Far[0].AbsLoc)
								 Reads[ReadIndex].UP_Far.clear();
						 }
						 else { // if (Reads[ReadIndex].MatchedD == Minus)
							 if (Reads[ReadIndex].UP_Far[0].AbsLoc < Reads[ReadIndex].ReadLength + Reads[ReadIndex].UP_Close[0].AbsLoc)
								 Reads[ReadIndex].UP_Far.clear();
						 }
					 }
					 else {
						 if (Reads[ReadIndex].UP_Far_backup.size()) {
							 if (Reads[ReadIndex].UP_Far_backup[Reads[ReadIndex].UP_Far_backup.size() - 1].LengthStr <Reads[ReadIndex].UP_Far[Reads[ReadIndex].UP_Far.size() - 1].LengthStr) {
								 Reads[ReadIndex].UP_Far_backup = Reads[ReadIndex].UP_Far;
								 Reads[ReadIndex].UP_Far.clear();
							 }
							 else Reads[ReadIndex].UP_Far.clear(); //CurrentChr
						 }
						 else {
							 Reads[ReadIndex].UP_Far_backup = Reads[ReadIndex].UP_Far;
							 Reads[ReadIndex].UP_Far.clear();
						 }
					 }
					 if (!Reads[ReadIndex].UP_Far.empty()) {
						 CountFarEnd++;
						 if (Reads[ReadIndex].MatchedD == Plus) CountFarEndPlus++;
						 else CountFarEndMinus++;
					 }
				 }
			 }

			 cout << RangeIndex << "\tNumber of reads with far end mapped: " << CountFarEnd << "\t"
			 //<< "\tTotal number of reads:" << Reads.size() << "\n"
			 //<< CountFarEnd * 100.0 / Reads.size() << " %\n"
			 << "Far+: " << CountFarEndPlus << "\tFar-: " << CountFarEndMinus << endl;
		 }

		 cout << "Searching breakpoints of inversions" << endl;
		 unsigned ReadsUsedForD = 0;
		 unsigned ReadsUsedForDI = 0;
		 for (short RangeIndex = 1; RangeIndex < MaxRangeIndex; RangeIndex++) {

			 CountFarEnd = 0;
			 CountFarEndMinus = 0;
			 CountFarEndPlus = 0;
			 for (unsigned ReadIndex = 0; ReadIndex < Reads.size(); ReadIndex++) {
				 if (Reads[ReadIndex].Used || Reads[ReadIndex].UP_Far.size()) {
					 continue;
				 }

				 //MinClose = short(log((double)(2 * Reads[ReadIndex].InsertSize + DSizeArray[RangeIndex]))/log(4.0) + 0.8) + 3;//atoi(argv[1]);
				 //MinFar_I = MinClose + 1;//atoi(argv[2]);
				 //cout << "For short insertion: " << MinClose << "\t" << MinFar_I << endl;
				 //MinFar_D = 8;//atoi(argv[3]);
				 //if (RangeIndex <= 5)
				   GetFarEnd_OtherStrand(CurrentChr, Reads[ReadIndex], RangeIndex);

				 if (!Reads[ReadIndex].UP_Far.empty()) {
					 if (Reads[ReadIndex].UP_Far[Reads[ReadIndex].UP_Far.size() - 1].LengthStr + Reads[ReadIndex].CloseEndLength < Reads[ReadIndex].ReadLength) {
						 if (Reads[ReadIndex].UP_Far_backup.size()) {
							 if (Reads[ReadIndex].UP_Far_backup[Reads[ReadIndex].UP_Far_backup.size() - 1].LengthStr <Reads[ReadIndex].UP_Far[Reads[ReadIndex].UP_Far.size() - 1].LengthStr) {
								 Reads[ReadIndex].UP_Far_backup = Reads[ReadIndex].UP_Far;
								 Reads[ReadIndex].UP_Far.clear();
							 }
							 else Reads[ReadIndex].UP_Far.clear();
						 }
						 else {
							 Reads[ReadIndex].UP_Far_backup = Reads[ReadIndex].UP_Far;
							 Reads[ReadIndex].UP_Far.clear();
						 }
					 }
					 else {
						 CountFarEnd++;
						 if (Reads[ReadIndex].MatchedD == Plus) CountFarEndPlus++;
						 else CountFarEndMinus++;
					 }
				 }
			 }
			 cout << RangeIndex << "\tNumber of reads with far end mapped: " << CountFarEnd << "\t"
			 //<< "\tTotal number of reads:" << Reads.size() << "\n"
			 //<< CountFarEnd * 100.0 / Reads.size() << " %\n"
			 << "Far+: " << CountFarEndPlus << "\tFar-: " << CountFarEndMinus << endl;
		 }
} // if (Analyze_TD_INV_LI_Others)

		 // compare backup with current value
		 cout << "revisit all breakpoints identified ...";
		 for (unsigned ReadIndex = 0; ReadIndex < Reads.size(); ReadIndex++) {
			 if (Reads[ReadIndex].UP_Far.empty()) {
				 if (!Reads[ReadIndex].UP_Far_backup.empty()) {
					 Reads[ReadIndex].UP_Far = Reads[ReadIndex].UP_Far_backup;
				 }
			 }
			 else if (!Reads[ReadIndex].UP_Far_backup.empty()) {
				 if (Reads[ReadIndex].UP_Far_backup[Reads[ReadIndex].UP_Far_backup.size() - 1].LengthStr >
					 Reads[ReadIndex].UP_Far[Reads[ReadIndex].UP_Far.size() - 1].LengthStr) {
					 Reads[ReadIndex].UP_Far = Reads[ReadIndex].UP_Far_backup;
				 }
			 }
		 }
		 cout << " done." << endl;


	    // ################### module 4: search variants and report #####################

		 //short MAX_MISMATCHES_Per_Read = 0;;
		 cout << "Searching deletion events ... " << endl;
		 for (unsigned ReadIndex = 0; ReadIndex < Reads.size(); ReadIndex++) {
			 if (Reads[ReadIndex].UP_Far.empty()) continue;
			 //MAX_MISMATCHES_Per_Read = (short)(Seq_Error_Rate * Reads[ReadIndex].ReadLength + 1);
			 //if (Reads[ReadIndex].UP_Far.size())
			 {
				 if (Reads[ReadIndex].MatchedD == Plus) { // MAX_SNP_ERROR
					 for (short MAX_SNP_ERROR_index = 0; MAX_SNP_ERROR_index <= Reads[ReadIndex].MAX_SNP_ERROR; MAX_SNP_ERROR_index++) {
					    for (int CloseIndex = 0; CloseIndex < Reads[ReadIndex].UP_Close.size(); CloseIndex++) {
							 if (Reads[ReadIndex].Used) break;
							 if (Reads[ReadIndex].UP_Close[CloseIndex].Mismatches > MAX_SNP_ERROR_index) continue;
							 for (int FarIndex = Reads[ReadIndex].UP_Far.size() - 1; FarIndex >= 0; FarIndex--) {
								 //cout << "+" << FarIndex << endl;
								 if (Reads[ReadIndex].Used) break;
								 if (Reads[ReadIndex].UP_Far[FarIndex].Mismatches > MAX_SNP_ERROR_index) continue;
								 if (Reads[ReadIndex].UP_Far[FarIndex].Mismatches + Reads[ReadIndex].UP_Close[CloseIndex].Mismatches > MAX_SNP_ERROR_index)
									 continue;
								 if (Reads[ReadIndex].UP_Far[FarIndex].Direction == Minus) {
									 //cout << "1" << endl;
									 //cout << "+\t" << Reads[ReadIndex].UP_Far[FarIndex].LengthStr
									 //	<< "\t" << Reads[ReadIndex].UP_Close[CloseIndex].LengthStr
									 //     << "\t" << Reads[ReadIndex].UP_Far[FarIndex].LengthStr + Reads[ReadIndex].UP_Close[CloseIndex].LengthStr << "\t" << Reads[ReadIndex].ReadLength << "\t";
									 //cout << Reads[ReadIndex].UP_Far[FarIndex].AbsLoc << "\t>\t" << Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc + 1 << endl;
									 if (Reads[ReadIndex].UP_Far[FarIndex].LengthStr + Reads[ReadIndex].UP_Close[CloseIndex].LengthStr == Reads[ReadIndex].ReadLength && Reads[ReadIndex].UP_Far[FarIndex].AbsLoc > Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc + 1) {
										 Reads[ReadIndex].Left = Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc - Reads[ReadIndex].UP_Close[CloseIndex].LengthStr + 1;
										 Reads[ReadIndex].Right = Reads[ReadIndex].UP_Far[FarIndex].AbsLoc + Reads[ReadIndex].UP_Far[FarIndex].LengthStr - 1;
										 Reads[ReadIndex].BP = Reads[ReadIndex].UP_Close[CloseIndex].LengthStr - 1;

										 Reads[ReadIndex].IndelSize =  (Reads[ReadIndex].Right - Reads[ReadIndex].Left) - Reads[ReadIndex].ReadLengthMinus;
										 Reads[ReadIndex].NT_str = "";
										 Reads[ReadIndex].NT_size = 0;
										 Reads[ReadIndex].InsertedStr = "";
										 Reads[ReadIndex].BPLeft = Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc - SpacerBeforeAfter;// - Fragments[LeftReads[Left_Index].MatchedSeqID].Start;
										 Reads[ReadIndex].BPRight = Reads[ReadIndex].UP_Far[FarIndex].AbsLoc - SpacerBeforeAfter;// - Fragments[LeftReads[Left_Index].MatchedSeqID].Start;
										 Reads[ReadIndex].score = Const_I + Const_S + LOG14 * Reads[ReadIndex].ReadLength + Const_Log_T;
										 //LeftReads[Left_Index].OK = true;
                                         unsigned RealBP_left= Reads[ReadIndex].BPLeft;
                                         unsigned RealBP_right = Reads[ReadIndex].BPRight;//, DIFF;
                                         if (Reads[ReadIndex].NT_str.size()) {
                                             GetRealStart4Insertion(CurrentChr, Reads[ReadIndex].NT_str, RealBP_left, RealBP_right);
                                         }
                                         else {
                                             GetRealStart4Deletion(CurrentChr, RealBP_left, RealBP_right);
                                         }
                                         short DIFF = Reads[ReadIndex].BPLeft - RealBP_left;
                                         DIFF = !((Reads[ReadIndex].BP - 1)<DIFF)?DIFF:(Reads[ReadIndex].BP - 1); // min(DIFF, currentRead.BP - 1);
                                         if (DIFF) {
                                             //std::cout << DIFF << std::endl;
                                             Reads[ReadIndex].BP -= DIFF;
                                             Reads[ReadIndex].BPLeft -= DIFF;
                                             Reads[ReadIndex].BPRight  -= DIFF;
                                         }
										 {
											 //if (LeftReads[Left_Index].IndelSize <= DSizeArray[RangeIndex] && LeftReads[Left_Index].IndelSize > DSizeArray[RangeIndex - 1])
											 {
												 Deletions[Reads[ReadIndex].BPLeft / BoxSize].push_back(ReadIndex);
												 Reads[ReadIndex].Used = true;
												 Count_D++;
												 Count_D_Plus++;
											 }
										 }
									 }
								 }
							 }
						 }
				    }
				 }
				 else if (Reads[ReadIndex].MatchedD == Minus) {
					 for (short MAX_SNP_ERROR_index = 0; MAX_SNP_ERROR_index <= Reads[ReadIndex].MAX_SNP_ERROR; MAX_SNP_ERROR_index++) {
					    for (int CloseIndex = Reads[ReadIndex].UP_Close.size() - 1; CloseIndex >= 0; CloseIndex--) {
							 if (Reads[ReadIndex].Used) break;
							 if (Reads[ReadIndex].UP_Close[CloseIndex].Mismatches > MAX_SNP_ERROR_index) continue;
							 for (int FarIndex = 0; FarIndex < Reads[ReadIndex].UP_Far.size(); FarIndex++) {
								 if (Reads[ReadIndex].Used) break;
								 if (Reads[ReadIndex].UP_Far[FarIndex].Mismatches > MAX_SNP_ERROR_index) continue;
								 if (Reads[ReadIndex].UP_Far[FarIndex].Mismatches + Reads[ReadIndex].UP_Close[CloseIndex].Mismatches > MAX_SNP_ERROR_index)
									 continue;
								 if (Reads[ReadIndex].UP_Far[FarIndex].Direction == Plus) {
									 if (Reads[ReadIndex].UP_Close[CloseIndex].LengthStr + Reads[ReadIndex].UP_Far[FarIndex].LengthStr == Reads[ReadIndex].ReadLength
										  && Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc > Reads[ReadIndex].UP_Far[FarIndex].AbsLoc + 1) {

										 Reads[ReadIndex].Left = Reads[ReadIndex].UP_Far[FarIndex].AbsLoc - Reads[ReadIndex].UP_Far[FarIndex].LengthStr + 1;
										 Reads[ReadIndex].Right = Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc + Reads[ReadIndex].UP_Close[CloseIndex].LengthStr - 1;
										 Reads[ReadIndex].BP = Reads[ReadIndex].UP_Far[FarIndex].LengthStr - 1;

										 Reads[ReadIndex].IndelSize = (Reads[ReadIndex].Right - Reads[ReadIndex].Left) - Reads[ReadIndex].ReadLengthMinus;
										 Reads[ReadIndex].NT_str = "";
										 Reads[ReadIndex].NT_size = 0;
										 Reads[ReadIndex].InsertedStr = "";
										 Reads[ReadIndex].BPLeft = Reads[ReadIndex].UP_Far[FarIndex].AbsLoc - SpacerBeforeAfter;// - Fragments[LeftReads[Left_Index].MatchedSeqID].Start;
										 Reads[ReadIndex].BPRight = Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc - SpacerBeforeAfter;// - Fragments[LeftReads[Left_Index].MatchedSeqID].Start;
										 Reads[ReadIndex].score = Const_I + Const_S + LOG14 * Reads[ReadIndex].ReadLength + Const_Log_T;
									 //LeftReads[Left_Index].OK = true;
                                         unsigned RealBP_left= Reads[ReadIndex].BPLeft;
                                         unsigned RealBP_right = Reads[ReadIndex].BPRight;//, DIFF;
                                         if (Reads[ReadIndex].NT_str.size()) {
                                             GetRealStart4Insertion(CurrentChr, Reads[ReadIndex].NT_str, RealBP_left, RealBP_right);
                                         }
                                         else {
                                             GetRealStart4Deletion(CurrentChr, RealBP_left, RealBP_right);
                                         }
                                         short DIFF = Reads[ReadIndex].BPLeft - RealBP_left;
                                         DIFF = !((Reads[ReadIndex].BP - 1)<DIFF)?DIFF:(Reads[ReadIndex].BP - 1); // min(DIFF, currentRead.BP - 1);
                                         if (DIFF) {
                                             //std::cout << DIFF << std::endl;
                                             Reads[ReadIndex].BP -= DIFF;
                                             Reads[ReadIndex].BPLeft -= DIFF;
                                             Reads[ReadIndex].BPRight  -= DIFF;
                                         }

										 {
											 //if (LeftReads[Left_Index].IndelSize <= DSizeArray[RangeIndex] && LeftReads[Left_Index].IndelSize > DSizeArray[RangeIndex - 1])
											 {
												 Deletions[Reads[ReadIndex].BPLeft / BoxSize].push_back(ReadIndex);
												 Reads[ReadIndex].Used = true;

												 Count_D++;
												 Count_D_Minus++;
												 //cout << "- " << Count_D << endl;
											 }
										 }
									 }
								 }
							 }
						 }
					 } // for (short MAX_SNP_ERROR_index = 0; MAX_SNP_ERROR_index <= Reads[ReadIndex].MAX_SNP_ERROR; MAX_SNP_ERROR_index++)
				 }
			}
		 }
			cout << "Total: " << Count_D << "\t+" << Count_D_Plus << "\t-" << Count_D_Minus << endl;
			ofstream DeletionOutf(DeletionOutputFilename.c_str());
			SortOutputD(NumBoxes, CurrentChr, Reads, Deletions, DeletionOutf);
			//DeletionOutf.close();
			for (unsigned int i = 0; i < NumBoxes; i++) Deletions[i].clear();


		 cout << "Searching deletion-insertions ... " << endl;
		 unsigned CloseIndex, FarIndex;
		 for (unsigned ReadIndex = 0; ReadIndex < Reads.size(); ReadIndex++) {
			 if (Reads[ReadIndex].Used || Reads[ReadIndex].UP_Far.empty()) continue;
			 //if (Reads[ReadIndex].UP_Far.size())
			 {
				 CloseIndex = Reads[ReadIndex].UP_Close.size() - 1;
				 FarIndex = Reads[ReadIndex].UP_Far.size() - 1;
				 if (Reads[ReadIndex].UP_Far[FarIndex].Mismatches + Reads[ReadIndex].UP_Close[CloseIndex].Mismatches > (short)(1 + Seq_Error_Rate * (Reads[ReadIndex].UP_Far[FarIndex].LengthStr + Reads[ReadIndex].UP_Close[CloseIndex].LengthStr)))
					 continue;
				 if (Reads[ReadIndex].MatchedD == Plus) {
					 //for (unsigned CloseIndex = 0; CloseIndex < Reads[ReadIndex].UP_Close.size(); CloseIndex++)
					 {
						 //for (unsigned FarIndex = 0; FarIndex < Reads[ReadIndex].UP_Far.size(); FarIndex++)
						 {

							 if (Reads[ReadIndex].UP_Far[FarIndex].Direction == Minus) {
								 if (Reads[ReadIndex].UP_Far[FarIndex].LengthStr + Reads[ReadIndex].UP_Close[CloseIndex].LengthStr < Reads[ReadIndex].ReadLength
									  && Reads[ReadIndex].UP_Far[FarIndex].LengthStr + Reads[ReadIndex].UP_Close[CloseIndex].LengthStr >= Min_Num_Matched_Bases && Reads[ReadIndex].UP_Far[FarIndex].AbsLoc > Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc + 1) {
									 //DeletionOK = true;
									 Reads[ReadIndex].Left = Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc - Reads[ReadIndex].UP_Close[CloseIndex].LengthStr + 1;
									 Reads[ReadIndex].Right = Reads[ReadIndex].UP_Far[FarIndex].AbsLoc + Reads[ReadIndex].UP_Far[FarIndex].LengthStr - 1;
									 Reads[ReadIndex].BP = Reads[ReadIndex].UP_Close[CloseIndex].LengthStr - 1;
									 Reads[ReadIndex].NT_size = Reads[ReadIndex].ReadLength - Reads[ReadIndex].UP_Far[FarIndex].LengthStr - Reads[ReadIndex].UP_Close[CloseIndex].LengthStr;

									 Reads[ReadIndex].NT_str = ReverseComplement(Reads[ReadIndex].UnmatchedSeq).substr(Reads[ReadIndex].BP + 1, Reads[ReadIndex].NT_size);
									 Reads[ReadIndex].InsertedStr = "";

									 Reads[ReadIndex].IndelSize =  (Reads[ReadIndex].Right - Reads[ReadIndex].Left) + Reads[ReadIndex].NT_size - Reads[ReadIndex].ReadLengthMinus;

									 Reads[ReadIndex].BPLeft = Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc - SpacerBeforeAfter;// - Fragments[LeftReads[Left_Index].MatchedSeqID].Start;
									 Reads[ReadIndex].BPRight = Reads[ReadIndex].UP_Far[FarIndex].AbsLoc - SpacerBeforeAfter;// - Fragments[LeftReads[Left_Index].MatchedSeqID].Start;
                                     unsigned RealBP_left= Reads[ReadIndex].BPLeft;
                                     unsigned RealBP_right = Reads[ReadIndex].BPRight;//, DIFF;
                                     if (Reads[ReadIndex].NT_str.size()) {
                                         GetRealStart4Insertion(CurrentChr, Reads[ReadIndex].NT_str, RealBP_left, RealBP_right);
                                     }
                                     else {
                                         GetRealStart4Deletion(CurrentChr, RealBP_left, RealBP_right);
                                     }
                                     short DIFF = Reads[ReadIndex].BPLeft - RealBP_left;
                                     DIFF = !((Reads[ReadIndex].BP - 1)<DIFF)?DIFF:(Reads[ReadIndex].BP - 1); // min(DIFF, currentRead.BP - 1);
                                     if (DIFF) {
                                         //std::cout << DIFF << std::endl;
                                         Reads[ReadIndex].BP -= DIFF;
                                         Reads[ReadIndex].BPLeft -= DIFF;
                                         Reads[ReadIndex].BPRight  -= DIFF;
                                     }

									 {
										 if (Reads[ReadIndex].IndelSize >= MIN_IndelSize_NT && Reads[ReadIndex].IndelSize > Reads[ReadIndex].NT_size)
										 {
											 DI[Reads[ReadIndex].BPLeft / BoxSize].push_back(ReadIndex);
											 Reads[ReadIndex].Used = true;
											 Count_DI++;
											 Count_DI_Plus++;
										 }
									 }

								 }

								 // ############
								 //}
							 }
						 }
					 }
				 }
				 else if (Reads[ReadIndex].MatchedD == Minus) {
					 //for (unsigned CloseIndex = 0; CloseIndex < Reads[ReadIndex].UP_Close.size(); CloseIndex++)
					 {
						 //for (unsigned FarIndex = 0; FarIndex < Reads[ReadIndex].UP_Far.size(); FarIndex++)
						 {
							 if (Reads[ReadIndex].UP_Far[FarIndex].Direction == Plus) {
								 if (Reads[ReadIndex].UP_Close[CloseIndex].LengthStr + Reads[ReadIndex].UP_Far[FarIndex].LengthStr < Reads[ReadIndex].ReadLength
									  && Reads[ReadIndex].UP_Close[CloseIndex].LengthStr + Reads[ReadIndex].UP_Far[FarIndex].LengthStr >= Min_Num_Matched_Bases && Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc > Reads[ReadIndex].UP_Far[FarIndex].AbsLoc + 1) {
									 //DeletionOK = true;
									 Reads[ReadIndex].Left = Reads[ReadIndex].UP_Far[FarIndex].AbsLoc - Reads[ReadIndex].UP_Far[FarIndex].LengthStr + 1;
									 Reads[ReadIndex].Right = Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc + Reads[ReadIndex].UP_Close[CloseIndex].LengthStr - 1;
									 Reads[ReadIndex].BP = Reads[ReadIndex].UP_Far[FarIndex].LengthStr - 1;
									 Reads[ReadIndex].NT_size = Reads[ReadIndex].ReadLength - Reads[ReadIndex].UP_Close[CloseIndex].LengthStr - Reads[ReadIndex].UP_Far[FarIndex].LengthStr;
									 Reads[ReadIndex].NT_str = Reads[ReadIndex].UnmatchedSeq.substr(Reads[ReadIndex].BP + 1, Reads[ReadIndex].NT_size);

									 Reads[ReadIndex].IndelSize =  (Reads[ReadIndex].Right - Reads[ReadIndex].Left) - Reads[ReadIndex].ReadLengthMinus + Reads[ReadIndex].NT_size;
									 //                                 LeftReads[Left_Index].NT_str = "";
									 Reads[ReadIndex].InsertedStr = "";
									 //cout << LeftReads[Left_Index].IndelSize << endl;
									 Reads[ReadIndex].BPLeft = Reads[ReadIndex].UP_Far[FarIndex].AbsLoc - SpacerBeforeAfter;// - Fragments[LeftReads[Left_Index].MatchedSeqID].Start;
									 Reads[ReadIndex].BPRight = Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc - SpacerBeforeAfter;// - Fragments[LeftReads[Left_Index].MatchedSeqID].Start;
									 //LeftReads[Left_Index].score = Const_I + Const_S + LOG14 * LeftReads[Left_Index].ReadLength + Const_Log_T;
									 //LeftReads[Left_Index].OK = true;
                                     unsigned RealBP_left= Reads[ReadIndex].BPLeft;
                                     unsigned RealBP_right = Reads[ReadIndex].BPRight;//, DIFF;
                                     if (Reads[ReadIndex].NT_str.size()) {
                                         GetRealStart4Insertion(CurrentChr, Reads[ReadIndex].NT_str, RealBP_left, RealBP_right);
                                     }
                                     else {
                                         GetRealStart4Deletion(CurrentChr, RealBP_left, RealBP_right);
                                     }
                                     short DIFF = Reads[ReadIndex].BPLeft - RealBP_left;
                                     DIFF = !((Reads[ReadIndex].BP - 1)<DIFF)?DIFF:(Reads[ReadIndex].BP - 1); // min(DIFF, currentRead.BP - 1);
                                     if (DIFF) {
                                         //std::cout << DIFF << std::endl;
                                         Reads[ReadIndex].BP -= DIFF;
                                         Reads[ReadIndex].BPLeft -= DIFF;
                                         Reads[ReadIndex].BPRight  -= DIFF;
                                     }

									 {
										 if (Reads[ReadIndex].IndelSize >= MIN_IndelSize_NT && Reads[ReadIndex].IndelSize > Reads[ReadIndex].NT_size)
										 {
											 DI[Reads[ReadIndex].BPLeft / BoxSize].push_back(ReadIndex);
											 Reads[ReadIndex].Used = true;
											 Count_DI++;
											 Count_DI_Minus++;
										 }
									 }

								 }

								 // ######################
							 }
						 }
					 }
				 }
			 }
		 }
		 cout << "Total: " << Count_DI << "\t+" << Count_DI_Plus << "\t-" << Count_DI_Minus << endl;
		 //ofstream DeletionInsertionOutf(DeletionInsertinOutputFilename.c_str());
		 SortOutputDI(NumBoxes, CurrentChr, Reads, DI, DeletionOutf);
		 //DeletionInsertionOutf.close();
		 DeletionOutf.close();
		 for (unsigned int i = 0; i < NumBoxes; i++) DI[i].clear();

if (Analyze_TD_INV_LI_Others) {


		 cout << "Searching tandem dupliation events ... " << endl;
		 for (unsigned ReadIndex = 0; ReadIndex < Reads.size(); ReadIndex++) {
			 if (Reads[ReadIndex].Used || Reads[ReadIndex].UP_Far.empty()) continue;
			 //short MAX_MISMATCHES_Per_Read = (short)(Seq_Error_Rate * Reads[ReadIndex].ReadLength + 1);
			 //if (Reads[ReadIndex].UP_Far.size())
			 {
				 if (Reads[ReadIndex].MatchedD == Plus) {
					 for (int CloseIndex = 0; CloseIndex < Reads[ReadIndex].UP_Close.size(); CloseIndex++) {
						 if (Reads[ReadIndex].Used) break;
						 for (int FarIndex = Reads[ReadIndex].UP_Far.size() - 1; FarIndex >= 0; FarIndex--) {
							 if (Reads[ReadIndex].Used) break;
							 if (Reads[ReadIndex].UP_Far[FarIndex].Mismatches + Reads[ReadIndex].UP_Close[CloseIndex].Mismatches > Reads[ReadIndex].MAX_SNP_ERROR)
								 continue;
							 if (Reads[ReadIndex].UP_Far[FarIndex].Direction == Minus) {
								 if (Reads[ReadIndex].UP_Far[FarIndex].LengthStr + Reads[ReadIndex].UP_Close[CloseIndex].LengthStr == Reads[ReadIndex].ReadLength && Reads[ReadIndex].UP_Far[FarIndex].AbsLoc + Reads[ReadIndex].ReadLength < Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc) {
									 Reads[ReadIndex].Right = Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc - Reads[ReadIndex].UP_Close[CloseIndex].LengthStr + 1;
									 Reads[ReadIndex].Left = Reads[ReadIndex].UP_Far[FarIndex].AbsLoc + Reads[ReadIndex].UP_Far[FarIndex].LengthStr - 1;
									 Reads[ReadIndex].BP = Reads[ReadIndex].UP_Close[CloseIndex].LengthStr - 1;

									 Reads[ReadIndex].IndelSize =  Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc - Reads[ReadIndex].UP_Far[FarIndex].AbsLoc + 1;
									 Reads[ReadIndex].NT_str = "";
									 Reads[ReadIndex].NT_size = 0;
									 Reads[ReadIndex].InsertedStr = "";
									 Reads[ReadIndex].BPRight = Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc - SpacerBeforeAfter;// - Fragments[LeftReads[Left_Index].MatchedSeqID].Start;
									 Reads[ReadIndex].BPLeft = Reads[ReadIndex].UP_Far[FarIndex].AbsLoc - SpacerBeforeAfter;// - Fragments[LeftReads[Left_Index].MatchedSeqID].Start;
									 Reads[ReadIndex].score = Const_I + Const_S + LOG14 * Reads[ReadIndex].ReadLength + Const_Log_T;
									 //LeftReads[Left_Index].OK = true;

									 {
										 //if (LeftReads[Left_Index].IndelSize <= DSizeArray[RangeIndex] && LeftReads[Left_Index].IndelSize > DSizeArray[RangeIndex - 1])
										 {
											 TD[Reads[ReadIndex].BPLeft / BoxSize].push_back(ReadIndex);
											 Reads[ReadIndex].Used = true;
											 Count_TD++;
											 Count_TD_Plus++;
										 }
									 }
								 }
							 }
						 }
					 }
				 }
				 else if (Reads[ReadIndex].MatchedD == Minus) {
					 for (int CloseIndex = Reads[ReadIndex].UP_Close.size() - 1; CloseIndex >= 0; CloseIndex--) {
						 if (Reads[ReadIndex].Used) break;
						 for (int FarIndex = 0; FarIndex < Reads[ReadIndex].UP_Far.size(); FarIndex++) {
							 if (Reads[ReadIndex].Used) break;
							 if (Reads[ReadIndex].UP_Far[FarIndex].Mismatches + Reads[ReadIndex].UP_Close[CloseIndex].Mismatches > Reads[ReadIndex].MAX_SNP_ERROR)
								 continue;
							 if (Reads[ReadIndex].UP_Far[FarIndex].Direction == Plus) {
								 if (Reads[ReadIndex].UP_Close[CloseIndex].LengthStr + Reads[ReadIndex].UP_Far[FarIndex].LengthStr == Reads[ReadIndex].ReadLength
									  && Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc + Reads[ReadIndex].ReadLength < Reads[ReadIndex].UP_Far[FarIndex].AbsLoc) {
									 Reads[ReadIndex].Right = Reads[ReadIndex].UP_Far[FarIndex].AbsLoc - Reads[ReadIndex].UP_Far[FarIndex].LengthStr + 1;
									 Reads[ReadIndex].Left = Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc + Reads[ReadIndex].UP_Close[CloseIndex].LengthStr - 1;
									 Reads[ReadIndex].BP = Reads[ReadIndex].UP_Far[FarIndex].LengthStr - 1;

									 Reads[ReadIndex].IndelSize = Reads[ReadIndex].UP_Far[FarIndex].AbsLoc - Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc + 1;;
									 Reads[ReadIndex].NT_str = "";
									 Reads[ReadIndex].NT_size = 0;
									 Reads[ReadIndex].InsertedStr = "";
									 Reads[ReadIndex].BPRight = Reads[ReadIndex].UP_Far[FarIndex].AbsLoc - SpacerBeforeAfter;// - Fragments[LeftReads[Left_Index].MatchedSeqID].Start;
									 Reads[ReadIndex].BPLeft = Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc - SpacerBeforeAfter;// - Fragments[LeftReads[Left_Index].MatchedSeqID].Start;
									 Reads[ReadIndex].score = Const_I + Const_S + LOG14 * Reads[ReadIndex].ReadLength + Const_Log_T;
									 //LeftReads[Left_Index].OK = true;

									 {
										 //if (LeftReads[Left_Index].IndelSize <= DSizeArray[RangeIndex] && LeftReads[Left_Index].IndelSize > DSizeArray[RangeIndex - 1])
										 {
											 TD[Reads[ReadIndex].BPLeft / BoxSize].push_back(ReadIndex);
											 Reads[ReadIndex].Used = true;

											 Count_TD++;
											 Count_TD_Minus++;
											 //cout << "- " << Count_D << endl;
										 }
									 }
								 }
							 }
						 }
					 }
				 }
			 }
		 }
		 cout << "Total: " << Count_TD << "\t+" << Count_TD_Plus << "\t-" << Count_TD_Minus << endl;
		 ofstream TDOutf(TDOutputFilename.c_str());
		 SortOutputTD(NumBoxes, CurrentChr, Reads, TD, TDOutf);
		 //TDOutf.close();
       for (unsigned int i = 0; i < NumBoxes; i++) TD[i].clear();

		 cout << "Searching tandem dupliation events with non-template sequence ... " << endl;
		 for (unsigned ReadIndex = 0; ReadIndex < Reads.size(); ReadIndex++) {
			 if (Reads[ReadIndex].Used || Reads[ReadIndex].UP_Far.empty()) continue;
			 //MAX_MISMATCHES_Per_Read = (short)(Seq_Error_Rate * Reads[ReadIndex].ReadLength + 1);
			 CloseIndex = Reads[ReadIndex].UP_Close.size() - 1;
			 FarIndex = Reads[ReadIndex].UP_Far.size() - 1;
			 if (Reads[ReadIndex].UP_Far[FarIndex].Mismatches + Reads[ReadIndex].UP_Close[CloseIndex].Mismatches > (short)(1 + Seq_Error_Rate * (Reads[ReadIndex].UP_Far[FarIndex].LengthStr + Reads[ReadIndex].UP_Close[CloseIndex].LengthStr)))
				 continue;
			 //if (Reads[ReadIndex].UP_Far.size())
			 {
				 if (Reads[ReadIndex].MatchedD == Plus) {
					 //for (int CloseIndex = 0; CloseIndex < Reads[ReadIndex].UP_Close.size(); CloseIndex++)
					 {
						 //if (Reads[ReadIndex].Used) break;
						 //for (int FarIndex = Reads[ReadIndex].UP_Far.size() - 1; FarIndex >= 0; FarIndex--)
						 {
							 //if (Reads[ReadIndex].Used) break;

							 if (Reads[ReadIndex].UP_Far[FarIndex].Direction == Minus) {
								 if (Reads[ReadIndex].UP_Far[FarIndex].LengthStr + Reads[ReadIndex].UP_Close[CloseIndex].LengthStr < Reads[ReadIndex].ReadLength && Reads[ReadIndex].UP_Far[FarIndex].AbsLoc + Reads[ReadIndex].ReadLength < Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc
									  && Reads[ReadIndex].UP_Far[FarIndex].LengthStr + Reads[ReadIndex].UP_Close[CloseIndex].LengthStr > Min_Num_Matched_Bases) {
									 Reads[ReadIndex].Right = Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc - Reads[ReadIndex].UP_Close[CloseIndex].LengthStr + 1;
									 Reads[ReadIndex].Left = Reads[ReadIndex].UP_Far[FarIndex].AbsLoc + Reads[ReadIndex].UP_Far[FarIndex].LengthStr - 1;
									 Reads[ReadIndex].BP = Reads[ReadIndex].UP_Close[CloseIndex].LengthStr - 1;

									 Reads[ReadIndex].IndelSize =  Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc - Reads[ReadIndex].UP_Far[FarIndex].AbsLoc + 1;
									 Reads[ReadIndex].NT_size = Reads[ReadIndex].ReadLength - Reads[ReadIndex].UP_Close[CloseIndex].LengthStr - Reads[ReadIndex].UP_Far[FarIndex].LengthStr;
									 Reads[ReadIndex].NT_str = ReverseComplement(Reads[ReadIndex].UnmatchedSeq).substr(Reads[ReadIndex].BP + 1, Reads[ReadIndex].NT_size);
									 Reads[ReadIndex].InsertedStr = "";
									 Reads[ReadIndex].BPRight = Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc - SpacerBeforeAfter;// - Fragments[LeftReads[Left_Index].MatchedSeqID].Start;
									 Reads[ReadIndex].BPLeft = Reads[ReadIndex].UP_Far[FarIndex].AbsLoc - SpacerBeforeAfter;// - Fragments[LeftReads[Left_Index].MatchedSeqID].Start;
									 Reads[ReadIndex].score = Const_I + Const_S + LOG14 * Reads[ReadIndex].ReadLength + Const_Log_T;
									 //LeftReads[Left_Index].OK = true;

									 {
										 //if (LeftReads[Left_Index].IndelSize <= DSizeArray[RangeIndex] && LeftReads[Left_Index].IndelSize > DSizeArray[RangeIndex - 1])
										 {
											 TD_NT[Reads[ReadIndex].BPLeft / BoxSize].push_back(ReadIndex);
											 Reads[ReadIndex].Used = true;
											 Count_TD_NT++;
											 Count_TD_NT_Plus++;
										 }
									 }
								 }
							 }
						 }
					 }
				 }
				 else if (Reads[ReadIndex].MatchedD == Minus) {
					 //for (int CloseIndex = Reads[ReadIndex].UP_Close.size() - 1; CloseIndex >= 0; CloseIndex--)
					 {
						 //if (Reads[ReadIndex].Used) break;
						 //for (int FarIndex = 0; FarIndex < Reads[ReadIndex].UP_Far.size(); FarIndex++)
						 {
							 if (Reads[ReadIndex].UP_Far[FarIndex].Direction == Plus) {
								 if (Reads[ReadIndex].UP_Close[CloseIndex].LengthStr + Reads[ReadIndex].UP_Far[FarIndex].LengthStr < Reads[ReadIndex].ReadLength
									  && Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc + Reads[ReadIndex].ReadLength < Reads[ReadIndex].UP_Far[FarIndex].AbsLoc
									  && Reads[ReadIndex].UP_Far[FarIndex].LengthStr + Reads[ReadIndex].UP_Close[CloseIndex].LengthStr > Min_Num_Matched_Bases) {

									 Reads[ReadIndex].Right = Reads[ReadIndex].UP_Far[FarIndex].AbsLoc - Reads[ReadIndex].UP_Far[FarIndex].LengthStr + 1;
									 Reads[ReadIndex].Left = Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc + Reads[ReadIndex].UP_Close[CloseIndex].LengthStr - 1;
									 Reads[ReadIndex].BP = Reads[ReadIndex].UP_Far[FarIndex].LengthStr - 1;

									 Reads[ReadIndex].IndelSize = Reads[ReadIndex].UP_Far[FarIndex].AbsLoc - Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc + 1;
									 Reads[ReadIndex].NT_size = Reads[ReadIndex].ReadLength - Reads[ReadIndex].UP_Close[CloseIndex].LengthStr - Reads[ReadIndex].UP_Far[FarIndex].LengthStr;
									 Reads[ReadIndex].NT_str = Reads[ReadIndex].UnmatchedSeq.substr(Reads[ReadIndex].BP + 1, Reads[ReadIndex].NT_size);
									 Reads[ReadIndex].InsertedStr = "";
									 Reads[ReadIndex].BPRight = Reads[ReadIndex].UP_Far[FarIndex].AbsLoc - SpacerBeforeAfter;// - Fragments[LeftReads[Left_Index].MatchedSeqID].Start;
									 Reads[ReadIndex].BPLeft = Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc - SpacerBeforeAfter;// - Fragments[LeftReads[Left_Index].MatchedSeqID].Start;
									 Reads[ReadIndex].score = Const_I + Const_S + LOG14 * Reads[ReadIndex].ReadLength + Const_Log_T;
									 //LeftReads[Left_Index].OK = true;

									 {
										 //if (LeftReads[Left_Index].IndelSize <= DSizeArray[RangeIndex] && LeftReads[Left_Index].IndelSize > DSizeArray[RangeIndex - 1])
										 {
											 TD_NT[Reads[ReadIndex].BPLeft / BoxSize].push_back(ReadIndex);
											 Reads[ReadIndex].Used = true;

											 Count_TD_NT++;
											 Count_TD_NT_Minus++;
											 //cout << "- " << Count_D << endl;
										 }
									 }
								 }
							 }
						 }
					 }
				 }
			 }
		 }
		 cout << "Total: " << Count_TD_NT << "\t+" << Count_TD_NT_Plus << "\t-" << Count_TD_NT_Minus << endl;
		 //ofstream TDOutf(TDOutputFilename.c_str());
		 SortOutputTD_NT(NumBoxes, CurrentChr, Reads, TD_NT, TDOutf);
		 TDOutf.close();
       for (unsigned int i = 0; i < NumBoxes; i++) TD_NT[i].clear();

		 cout << "Searching inversions ... " << endl;
		 for (unsigned ReadIndex = 0; ReadIndex < Reads.size(); ReadIndex++) {
			 if (Reads[ReadIndex].Used || Reads[ReadIndex].UP_Far.empty()) continue;
			 if (Reads[ReadIndex].UP_Close[0].Strand != Reads[ReadIndex].UP_Far[0].Strand
				  &&
				  Reads[ReadIndex].UP_Close[0].Direction == Reads[ReadIndex].UP_Far[0].Direction) {
				 if (Reads[ReadIndex].MatchedD == Plus) {
					 for (int CloseIndex = Reads[ReadIndex].UP_Close.size() - 1; CloseIndex >= 0; CloseIndex--) {
						 if (Reads[ReadIndex].Used) break;
						 for (int FarIndex = Reads[ReadIndex].UP_Far.size() - 1; FarIndex >= 0; FarIndex--) {
							 if (Reads[ReadIndex].Used) break;
							 if (Reads[ReadIndex].UP_Far[FarIndex].Mismatches + Reads[ReadIndex].UP_Close[CloseIndex].Mismatches > Reads[ReadIndex].MAX_SNP_ERROR)
								 continue;
							 if (Reads[ReadIndex].UP_Far[FarIndex].Direction == Plus) {
								 if (Reads[ReadIndex].UP_Far[FarIndex].LengthStr + Reads[ReadIndex].UP_Close[CloseIndex].LengthStr == Reads[ReadIndex].ReadLength && Reads[ReadIndex].UP_Far[FarIndex].AbsLoc > Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc + MIN_IndelSize_Inversion) {
									 //cout << "+\t" << Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc
									 //<< "\t" << Reads[ReadIndex].UP_Far[FarIndex].AbsLoc << endl;
									 Reads[ReadIndex].Left = (Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc + 1) - Reads[ReadIndex].UP_Close[CloseIndex].LengthStr;
									 Reads[ReadIndex].Right = Reads[ReadIndex].UP_Far[FarIndex].AbsLoc - Reads[ReadIndex].UP_Far[FarIndex].LengthStr + Reads[ReadIndex].ReadLength;
									 Reads[ReadIndex].BP = Reads[ReadIndex].UP_Close[CloseIndex].LengthStr - 1;

									 Reads[ReadIndex].IndelSize =  Reads[ReadIndex].UP_Far[FarIndex].AbsLoc - Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc;
									 Reads[ReadIndex].NT_str = "";
									 Reads[ReadIndex].NT_size = 0;
									 Reads[ReadIndex].InsertedStr = "";
									 Reads[ReadIndex].BPLeft = Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc + 1 - SpacerBeforeAfter;// - Fragments[LeftReads[Left_Index].MatchedSeqID].Start;
									 Reads[ReadIndex].BPRight = Reads[ReadIndex].UP_Far[FarIndex].AbsLoc - SpacerBeforeAfter;// - Fragments[LeftReads[Left_Index].MatchedSeqID].Start;
									 Reads[ReadIndex].score = Const_I + Const_S + LOG14 * Reads[ReadIndex].ReadLength + Const_Log_T;
									 //LeftReads[Left_Index].OK = true;

									 {
										 //if (LeftReads[Left_Index].IndelSize <= DSizeArray[RangeIndex] && LeftReads[Left_Index].IndelSize > DSizeArray[RangeIndex - 1])
										 {
											 Inv[Reads[ReadIndex].BPLeft / BoxSize].push_back(ReadIndex);
											 Reads[ReadIndex].Used = true;
											 Count_Inv++;
											 Count_Inv_Plus++;
										 }
									 }
								 }
								 // anchor inside reversed block.
								 if (Reads[ReadIndex].UP_Far[FarIndex].LengthStr + Reads[ReadIndex].UP_Close[CloseIndex].LengthStr == Reads[ReadIndex].ReadLength && Reads[ReadIndex].UP_Far[FarIndex].AbsLoc + MIN_IndelSize_Inversion < Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc
									  ) {
									 //cout << "+\t" << Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc
									 //<< "\t" << Reads[ReadIndex].UP_Far[FarIndex].AbsLoc << endl;
									 Reads[ReadIndex].Right = Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc - Reads[ReadIndex].UP_Close[CloseIndex].LengthStr + Reads[ReadIndex].ReadLength;
									 Reads[ReadIndex].Left = Reads[ReadIndex].UP_Far[FarIndex].AbsLoc - Reads[ReadIndex].UP_Far[FarIndex].LengthStr + 1;
									 Reads[ReadIndex].BP = Reads[ReadIndex].UP_Far[FarIndex].LengthStr - 1;

									 Reads[ReadIndex].IndelSize = Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc - Reads[ReadIndex].UP_Far[FarIndex].AbsLoc;
									 Reads[ReadIndex].NT_str = "";
									 Reads[ReadIndex].NT_size = 0;
									 Reads[ReadIndex].InsertedStr = "";
									 Reads[ReadIndex].BPRight = Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc - SpacerBeforeAfter;// - Fragments[LeftReads[Left_Index].MatchedSeqID].Start;
									 Reads[ReadIndex].BPLeft = (Reads[ReadIndex].UP_Far[FarIndex].AbsLoc + 1) - SpacerBeforeAfter;// - Fragments[LeftReads[Left_Index].MatchedSeqID].Start;
									 Reads[ReadIndex].score = Const_I + Const_S + LOG14 * Reads[ReadIndex].ReadLength + Const_Log_T;
									 //LeftReads[Left_Index].OK = true;

									 {
										 //if (LeftReads[Left_Index].IndelSize <= DSizeArray[RangeIndex] && LeftReads[Left_Index].IndelSize > DSizeArray[RangeIndex - 1])
										 {
											 Inv[Reads[ReadIndex].BPLeft / BoxSize].push_back(ReadIndex);
											 Reads[ReadIndex].Used = true;
											 Count_Inv++;
											 Count_Inv_Plus++;
										 }
									 }
								 }
							 }
						 }
					 }
				 }
				 else if (Reads[ReadIndex].MatchedD == Minus) {
					 for (int CloseIndex = Reads[ReadIndex].UP_Close.size() - 1; CloseIndex >= 0; CloseIndex--) {
						 if (Reads[ReadIndex].Used) break;
						 for (int FarIndex = Reads[ReadIndex].UP_Far.size() - 1; FarIndex >= 0; FarIndex--) {
							 if (Reads[ReadIndex].Used) break;
							 if (Reads[ReadIndex].UP_Far[FarIndex].Mismatches + Reads[ReadIndex].UP_Close[CloseIndex].Mismatches > Reads[ReadIndex].MAX_SNP_ERROR)
								 continue;
							 if (Reads[ReadIndex].UP_Far[FarIndex].Direction == Minus) {
								 // ######################
								 //cout << "-\t" << Reads[ReadIndex].UP_Close[CloseIndex].LengthStr + Reads[ReadIndex].UP_Far[FarIndex].LengthStr << "\t" << Reads[ReadIndex].ReadLength << "\t";
								 //cout << Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc << "\t>\t" << Reads[ReadIndex].UP_Far[FarIndex].AbsLoc + 1 << endl;
								 // anchor outside reversed block.
								 if (Reads[ReadIndex].UP_Close[CloseIndex].LengthStr + Reads[ReadIndex].UP_Far[FarIndex].LengthStr == Reads[ReadIndex].ReadLength
									  && Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc > Reads[ReadIndex].UP_Far[FarIndex].AbsLoc + MIN_IndelSize_Inversion) {
									 //cout << "-\t" << Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc
									 //<< "\t" << Reads[ReadIndex].UP_Far[FarIndex].AbsLoc << endl;
									 Reads[ReadIndex].Left = Reads[ReadIndex].UP_Far[FarIndex].AbsLoc + Reads[ReadIndex].UP_Far[FarIndex].LengthStr - Reads[ReadIndex].ReadLength;
									 Reads[ReadIndex].Right = Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc + Reads[ReadIndex].UP_Close[CloseIndex].LengthStr - 1;
									 Reads[ReadIndex].BP = Reads[ReadIndex].UP_Far[FarIndex].LengthStr - 1;

									 Reads[ReadIndex].IndelSize = Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc - Reads[ReadIndex].UP_Far[FarIndex].AbsLoc;
									 Reads[ReadIndex].NT_str = "";
									 Reads[ReadIndex].NT_size = 0;
									 Reads[ReadIndex].InsertedStr = "";
									 Reads[ReadIndex].BPLeft = Reads[ReadIndex].UP_Far[FarIndex].AbsLoc - SpacerBeforeAfter;// - Fragments[LeftReads[Left_Index].MatchedSeqID].Start;
									 //cout <<  "far\t" << Reads[ReadIndex].UP_Far[FarIndex].AbsLoc << "\tspacer\t" << SpacerBeforeAfter << endl;
									 //cout <<  Reads[ReadIndex].UP_Far[FarIndex].AbsLoc - SpacerBeforeAfter - 1 << endl;
									 Reads[ReadIndex].BPRight = Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc - 1 - SpacerBeforeAfter;// - Fragments[LeftReads[Left_Index].MatchedSeqID].Start;
									 Reads[ReadIndex].score = Const_I + Const_S + LOG14 * Reads[ReadIndex].ReadLength + Const_Log_T;
									 //LeftReads[Left_Index].OK = true;

									 {
										 //if (LeftReads[Left_Index].IndelSize <= DSizeArray[RangeIndex] && LeftReads[Left_Index].IndelSize > DSizeArray[RangeIndex - 1])
										 {
											 //cout << "Inv\t" << Reads[ReadIndex].BPLeft << "\t" << BoxSize << endl;
											 Inv[Reads[ReadIndex].BPLeft / BoxSize].push_back(ReadIndex);
											 Reads[ReadIndex].Used = true;

											 Count_Inv++;
											 Count_Inv_Minus++;

										 }
									 }
								 }
								 // anchor inside reversed block.
								 if (Reads[ReadIndex].UP_Close[CloseIndex].LengthStr + Reads[ReadIndex].UP_Far[FarIndex].LengthStr == Reads[ReadIndex].ReadLength
									  && Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc + MIN_IndelSize_Inversion < Reads[ReadIndex].UP_Far[FarIndex].AbsLoc
									  ) {
									 Reads[ReadIndex].Right = Reads[ReadIndex].UP_Far[FarIndex].AbsLoc + Reads[ReadIndex].UP_Far[FarIndex].LengthStr -  1;
									 Reads[ReadIndex].Left = Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc + Reads[ReadIndex].UP_Close[CloseIndex].LengthStr - Reads[ReadIndex].ReadLength;
									 Reads[ReadIndex].BP = Reads[ReadIndex].UP_Close[CloseIndex].LengthStr - 1;

									 Reads[ReadIndex].IndelSize = Reads[ReadIndex].UP_Far[FarIndex].AbsLoc - Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc;
									 Reads[ReadIndex].NT_str = "";
									 Reads[ReadIndex].NT_size = 0;
									 Reads[ReadIndex].InsertedStr = "";
									 Reads[ReadIndex].BPLeft = Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc - SpacerBeforeAfter;// - Fragments[LeftReads[Left_Index].MatchedSeqID].Start;
									 //cout <<  "far\t" << Reads[ReadIndex].UP_Far[FarIndex].AbsLoc << "\tspacer\t" << SpacerBeforeAfter << endl;
									 //cout <<  Reads[ReadIndex].UP_Far[FarIndex].AbsLoc - SpacerBeforeAfter - 1 << endl;
									 Reads[ReadIndex].BPRight = Reads[ReadIndex].UP_Far[FarIndex].AbsLoc - 1 - SpacerBeforeAfter;// - Fragments[LeftReads[Left_Index].MatchedSeqID].Start;
									 Reads[ReadIndex].score = Const_I + Const_S + LOG14 * Reads[ReadIndex].ReadLength + Const_Log_T;
									 //LeftReads[Left_Index].OK = true;

									 {
										 //if (LeftReads[Left_Index].IndelSize <= DSizeArray[RangeIndex] && LeftReads[Left_Index].IndelSize > DSizeArray[RangeIndex - 1])
										 {
											 //cout << "Inv\t" << Reads[ReadIndex].BPLeft << "\t" << BoxSize << endl;
											 Inv[Reads[ReadIndex].BPLeft / BoxSize].push_back(ReadIndex);
											 Reads[ReadIndex].Used = true;

											 Count_Inv++;
											 Count_Inv_Minus++;
											 /*
											 cout << "4" << endl;
											 cout << Cap2Low(ReverseComplement(CurrentChr.substr(Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc, ReportLength)));
											 cout << CurrentChr.substr(Reads[ReadIndex].UP_Far[FarIndex].AbsLoc, ReportLength) << endl;
											 for (int i = 0; i < ReportLength - Reads[ReadIndex].BP - 1; i++) cout << " ";
											 cout << ReverseComplement( Reads[ReadIndex].UnmatchedSeq ) << endl;

											 DeletionInsertionOutf << Reads[ReadIndex].Name << "\n"
											 << Reads[ReadIndex].UnmatchedSeq << "\n"
											 << Reads[ReadIndex].MatchedD << "\t"
											 << Reads[ReadIndex].FragName << "\t"
											 << Reads[ReadIndex].MatchedRelPos << "\t"
											 << Reads[ReadIndex].MS << "\t"
											 << Reads[ReadIndex].InsertSize << "\t"
											 << Reads[ReadIndex].Tag << endl;
											 //cout << "- " << Count_D << endl;
											  */
										 }
									 }
								 }
							 }
						 }
					 }
				 }
			 }
		 }
		 cout << "Total: " << Count_Inv << "\t+" << Count_Inv_Plus << "\t-" << Count_Inv_Minus << endl;
		 ofstream InversionOutf(InversionOutputFilename.c_str());
		 SortOutputInv(NumBoxes, CurrentChr, Reads, Inv, InversionOutf);
		 //InversionOutf.close();
       for (unsigned int i = 0; i < NumBoxes; i++) Inv[i].clear();

		 cout << "Searching inversions with non-template sequence ... " << endl;
		 for (unsigned ReadIndex = 0; ReadIndex < Reads.size(); ReadIndex++) {
			 if (Reads[ReadIndex].Used || Reads[ReadIndex].UP_Far.empty()) continue;
			 CloseIndex = Reads[ReadIndex].UP_Close.size() - 1;
			 FarIndex = Reads[ReadIndex].UP_Far.size() - 1;
			 if (Reads[ReadIndex].UP_Far[FarIndex].Mismatches + Reads[ReadIndex].UP_Close[CloseIndex].Mismatches > (short)(1 + Seq_Error_Rate * (Reads[ReadIndex].UP_Far[FarIndex].LengthStr + Reads[ReadIndex].UP_Close[CloseIndex].LengthStr)))
				 continue;
			 if (Reads[ReadIndex].UP_Close[0].Strand != Reads[ReadIndex].UP_Far[0].Strand
				  &&
				  Reads[ReadIndex].UP_Close[0].Direction == Reads[ReadIndex].UP_Far[0].Direction) {
				 if (Reads[ReadIndex].MatchedD == Plus) {
					 {
						 {
							 //if (Reads[ReadIndex].Used) break;

							 if (Reads[ReadIndex].UP_Far[FarIndex].Direction == Plus) {
								 if (Reads[ReadIndex].UP_Far[FarIndex].LengthStr + Reads[ReadIndex].UP_Close[CloseIndex].LengthStr < Reads[ReadIndex].ReadLength
									  &&
									  Reads[ReadIndex].UP_Far[FarIndex].AbsLoc > Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc + MIN_IndelSize_Inversion
									  && Reads[ReadIndex].UP_Far[FarIndex].LengthStr + Reads[ReadIndex].UP_Close[CloseIndex].LengthStr >= Min_Num_Matched_Bases) {
									 //cout << "+\t" << Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc
									 //<< "\t" << Reads[ReadIndex].UP_Far[FarIndex].AbsLoc << endl;
									 Reads[ReadIndex].Left = (Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc + 1) - Reads[ReadIndex].UP_Close[CloseIndex].LengthStr;
									 Reads[ReadIndex].Right = Reads[ReadIndex].UP_Far[FarIndex].AbsLoc - Reads[ReadIndex].UP_Far[FarIndex].LengthStr + Reads[ReadIndex].ReadLength;
									 Reads[ReadIndex].BP = Reads[ReadIndex].UP_Close[CloseIndex].LengthStr - 1;

									 Reads[ReadIndex].IndelSize =  Reads[ReadIndex].UP_Far[FarIndex].AbsLoc - Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc;

									 Reads[ReadIndex].NT_size = Reads[ReadIndex].ReadLength - Reads[ReadIndex].UP_Far[FarIndex].LengthStr - Reads[ReadIndex].UP_Close[CloseIndex].LengthStr; // NT_2str
									 //cout << "Po " << Reads[ReadIndex].NT_size << "\t" <<  Reads[ReadIndex].ReadLength << "\t" << Reads[ReadIndex].UP_Close[CloseIndex].LengthStr << "\t" << Reads[ReadIndex].UP_Far[FarIndex].LengthStr << endl;
									 Reads[ReadIndex].NT_str = ReverseComplement(Reads[ReadIndex].UnmatchedSeq).substr(Reads[ReadIndex].BP + 1, Reads[ReadIndex].NT_size);
									 Reads[ReadIndex].InsertedStr = "";
									 Reads[ReadIndex].BPLeft = Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc + 1 - SpacerBeforeAfter;// - Fragments[LeftReads[Left_Index].MatchedSeqID].Start;
									 Reads[ReadIndex].BPRight = Reads[ReadIndex].UP_Far[FarIndex].AbsLoc - SpacerBeforeAfter;// - Fragments[LeftReads[Left_Index].MatchedSeqID].Start;
									 Reads[ReadIndex].score = Const_I + Const_S + LOG14 * Reads[ReadIndex].ReadLength + Const_Log_T;
									 //LeftReads[Left_Index].OK = true;

									 {
										 //if (LeftReads[Left_Index].IndelSize <= DSizeArray[RangeIndex] && LeftReads[Left_Index].IndelSize > DSizeArray[RangeIndex - 1])
										 {
											 Inv_NT[Reads[ReadIndex].BPLeft / BoxSize].push_back(ReadIndex);
											 Reads[ReadIndex].Used = true;
											 Count_Inv_NT++;
											 Count_Inv_NT_Plus++;
										 }
									 }
								 }
								 // anchor inside reversed block.
								 if (Reads[ReadIndex].UP_Far[FarIndex].LengthStr + Reads[ReadIndex].UP_Close[CloseIndex].LengthStr < Reads[ReadIndex].ReadLength && Reads[ReadIndex].UP_Far[FarIndex].AbsLoc + MIN_IndelSize_Inversion < Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc
									  && Reads[ReadIndex].UP_Far[FarIndex].LengthStr + Reads[ReadIndex].UP_Close[CloseIndex].LengthStr >= Min_Num_Matched_Bases
									  ) {
									 //cout << "+\t" << Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc
									 //<< "\t" << Reads[ReadIndex].UP_Far[FarIndex].AbsLoc << endl;
									 Reads[ReadIndex].Right = Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc - Reads[ReadIndex].UP_Close[CloseIndex].LengthStr + Reads[ReadIndex].ReadLength;
									 Reads[ReadIndex].Left = Reads[ReadIndex].UP_Far[FarIndex].AbsLoc - Reads[ReadIndex].UP_Far[FarIndex].LengthStr + 1;
									 Reads[ReadIndex].BP = Reads[ReadIndex].UP_Far[FarIndex].LengthStr - 1;

									 Reads[ReadIndex].IndelSize = Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc - Reads[ReadIndex].UP_Far[FarIndex].AbsLoc;

									 Reads[ReadIndex].NT_size = Reads[ReadIndex].ReadLength - Reads[ReadIndex].UP_Far[FarIndex].LengthStr - Reads[ReadIndex].UP_Close[CloseIndex].LengthStr;
									 //cout << "Pi " << Reads[ReadIndex].NT_size << "\t" <<  Reads[ReadIndex].ReadLength << "\t" << Reads[ReadIndex].UP_Close[CloseIndex].LengthStr << "\t" << Reads[ReadIndex].UP_Far[FarIndex].LengthStr << endl;
									 Reads[ReadIndex].NT_str = Reads[ReadIndex].UnmatchedSeq.substr(Reads[ReadIndex].BP + 1, Reads[ReadIndex].NT_size);
									 Reads[ReadIndex].InsertedStr = "";
									 Reads[ReadIndex].BPRight = Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc - SpacerBeforeAfter;// - Fragments[LeftReads[Left_Index].MatchedSeqID].Start;
									 Reads[ReadIndex].BPLeft = (Reads[ReadIndex].UP_Far[FarIndex].AbsLoc + 1) - SpacerBeforeAfter;// - Fragments[LeftReads[Left_Index].MatchedSeqID].Start;
									 Reads[ReadIndex].score = Const_I + Const_S + LOG14 * Reads[ReadIndex].ReadLength + Const_Log_T;
									 //LeftReads[Left_Index].OK = true;

									 {
										 //if (LeftReads[Left_Index].IndelSize <= DSizeArray[RangeIndex] && LeftReads[Left_Index].IndelSize > DSizeArray[RangeIndex - 1])
										 {
											 Inv_NT[Reads[ReadIndex].BPLeft / BoxSize].push_back(ReadIndex);
											 Reads[ReadIndex].Used = true;
											 Count_Inv_NT++;
											 Count_Inv_NT_Plus++;
										 }
									 }
								 }
							 }
						 }
					 }
				 }
				 else if (Reads[ReadIndex].MatchedD == Minus) {
					 //for (int CloseIndex = Reads[ReadIndex].UP_Close.size() - 1; CloseIndex >= 0; CloseIndex--)
					 {
						 //if (Reads[ReadIndex].Used) break;
						 //for (int FarIndex = Reads[ReadIndex].UP_Far.size() - 1; FarIndex >= 0; FarIndex--)
						 {
							 //if (Reads[ReadIndex].Used) break;
							 if (Reads[ReadIndex].UP_Far[FarIndex].Direction == Minus) {
								 // ######################
								 //cout << "-\t" << Reads[ReadIndex].UP_Close[CloseIndex].LengthStr + Reads[ReadIndex].UP_Far[FarIndex].LengthStr << "\t" << Reads[ReadIndex].ReadLength << "\t";
								 //cout << Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc << "\t>\t" << Reads[ReadIndex].UP_Far[FarIndex].AbsLoc + 1 << endl;
								 // anchor outside reversed block.
								 if (Reads[ReadIndex].UP_Close[CloseIndex].LengthStr + Reads[ReadIndex].UP_Far[FarIndex].LengthStr < Reads[ReadIndex].ReadLength
									  && Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc > Reads[ReadIndex].UP_Far[FarIndex].AbsLoc + MIN_IndelSize_Inversion
									  && Reads[ReadIndex].UP_Far[FarIndex].LengthStr + Reads[ReadIndex].UP_Close[CloseIndex].LengthStr >= Min_Num_Matched_Bases) {
									 //cout << "-\t" << Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc
									 //<< "\t" << Reads[ReadIndex].UP_Far[FarIndex].AbsLoc << endl;
									 Reads[ReadIndex].Left = Reads[ReadIndex].UP_Far[FarIndex].AbsLoc + Reads[ReadIndex].UP_Far[FarIndex].LengthStr - Reads[ReadIndex].ReadLength;
									 Reads[ReadIndex].Right = Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc + Reads[ReadIndex].UP_Close[CloseIndex].LengthStr - 1;
									 Reads[ReadIndex].BP = Reads[ReadIndex].UP_Far[FarIndex].LengthStr - 1;

									 Reads[ReadIndex].IndelSize = Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc - Reads[ReadIndex].UP_Far[FarIndex].AbsLoc;

									 Reads[ReadIndex].NT_size = Reads[ReadIndex].ReadLength - Reads[ReadIndex].UP_Far[FarIndex].LengthStr - Reads[ReadIndex].UP_Close[CloseIndex].LengthStr;
									 //cout << "Mo " << Reads[ReadIndex].NT_2size << "\t" <<  Reads[ReadIndex].ReadLength << "\t" << Reads[ReadIndex].UP_Close[CloseIndex].LengthStr << "\t" << Reads[ReadIndex].UP_Far[FarIndex].LengthStr << endl;
									 Reads[ReadIndex].NT_str = Reads[ReadIndex].UnmatchedSeq.substr(Reads[ReadIndex].BP + 1, Reads[ReadIndex].NT_size);
									 Reads[ReadIndex].InsertedStr = "";
									 Reads[ReadIndex].BPLeft = Reads[ReadIndex].UP_Far[FarIndex].AbsLoc - SpacerBeforeAfter;// - Fragments[LeftReads[Left_Index].MatchedSeqID].Start;
									 //cout <<  "far\t" << Reads[ReadIndex].UP_Far[FarIndex].AbsLoc << "\tspacer\t" << SpacerBeforeAfter << endl;
									 //cout <<  Reads[ReadIndex].UP_Far[FarIndex].AbsLoc - SpacerBeforeAfter - 1 << endl;
									 Reads[ReadIndex].BPRight = Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc - 1 - SpacerBeforeAfter;// - Fragments[LeftReads[Left_Index].MatchedSeqID].Start;
									 Reads[ReadIndex].score = Const_I + Const_S + LOG14 * Reads[ReadIndex].ReadLength + Const_Log_T;
									 //LeftReads[Left_Index].OK = true;

									 {
										 //if (LeftReads[Left_Index].IndelSize <= DSizeArray[RangeIndex] && LeftReads[Left_Index].IndelSize > DSizeArray[RangeIndex - 1])
										 {
											 //cout << "Inv\t" << Reads[ReadIndex].BPLeft << "\t" << BoxSize << endl;
											 Inv_NT[Reads[ReadIndex].BPLeft / BoxSize].push_back(ReadIndex);
											 Reads[ReadIndex].Used = true;

											 Count_Inv_NT++;
											 Count_Inv_NT_Minus++;

										 }
									 }
								 }
								 // anchor inside reversed block.
								 if (Reads[ReadIndex].UP_Close[CloseIndex].LengthStr + Reads[ReadIndex].UP_Far[FarIndex].LengthStr < Reads[ReadIndex].ReadLength
									  && Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc + MIN_IndelSize_Inversion < Reads[ReadIndex].UP_Far[FarIndex].AbsLoc
									  && Reads[ReadIndex].UP_Far[FarIndex].LengthStr + Reads[ReadIndex].UP_Close[CloseIndex].LengthStr >= Min_Num_Matched_Bases
									  ) {
									 Reads[ReadIndex].Right = Reads[ReadIndex].UP_Far[FarIndex].AbsLoc + Reads[ReadIndex].UP_Far[FarIndex].LengthStr -  1;
									 Reads[ReadIndex].Left = Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc + Reads[ReadIndex].UP_Close[CloseIndex].LengthStr - Reads[ReadIndex].ReadLength;
									 Reads[ReadIndex].BP = Reads[ReadIndex].UP_Close[CloseIndex].LengthStr - 1;

									 Reads[ReadIndex].IndelSize = Reads[ReadIndex].UP_Far[FarIndex].AbsLoc - Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc;

									 Reads[ReadIndex].NT_size = Reads[ReadIndex].ReadLength - Reads[ReadIndex].UP_Far[FarIndex].LengthStr - Reads[ReadIndex].UP_Close[CloseIndex].LengthStr;
									 //cout << "Mi " << Reads[ReadIndex].NT_2size << "\t" <<  Reads[ReadIndex].ReadLength << "\t" << Reads[ReadIndex].UP_Close[CloseIndex].LengthStr << "\t" << Reads[ReadIndex].UP_Far[FarIndex].LengthStr << endl;
									 Reads[ReadIndex].NT_str = ReverseComplement(Reads[ReadIndex].UnmatchedSeq).substr(Reads[ReadIndex].BP + 1, Reads[ReadIndex].NT_size);
									 Reads[ReadIndex].InsertedStr = "";
									 Reads[ReadIndex].BPLeft = Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc - SpacerBeforeAfter;// - Fragments[LeftReads[Left_Index].MatchedSeqID].Start;
									 //cout <<  "far\t" << Reads[ReadIndex].UP_Far[FarIndex].AbsLoc << "\tspacer\t" << SpacerBeforeAfter << endl;
									 //cout <<  Reads[ReadIndex].UP_Far[FarIndex].AbsLoc - SpacerBeforeAfter - 1 << endl;
									 Reads[ReadIndex].BPRight = Reads[ReadIndex].UP_Far[FarIndex].AbsLoc - 1 - SpacerBeforeAfter;// - Fragments[LeftReads[Left_Index].MatchedSeqID].Start;
									 Reads[ReadIndex].score = Const_I + Const_S + LOG14 * Reads[ReadIndex].ReadLength + Const_Log_T;
									 //LeftReads[Left_Index].OK = true;

									 {
										 //if (LeftReads[Left_Index].IndelSize <= DSizeArray[RangeIndex] && LeftReads[Left_Index].IndelSize > DSizeArray[RangeIndex - 1])
										 {
											 //cout << "Inv\t" << Reads[ReadIndex].BPLeft << "\t" << BoxSize << endl;
											 Inv_NT[Reads[ReadIndex].BPLeft / BoxSize].push_back(ReadIndex);
											 Reads[ReadIndex].Used = true;

											 Count_Inv_NT++;
											 Count_Inv_NT_Minus++;
										 }
									 }
								 }
							 }
						 }
					 }
				 }
			 }
		 }
		 cout << "Total: " << Count_Inv_NT << "\t+" << Count_Inv_NT_Plus << "\t-" << Count_Inv_NT_Minus << endl;
		 //ofstream InversionOutf(InversionOutputFilename.c_str());
		 SortOutputInv_NT(NumBoxes, CurrentChr, Reads, Inv_NT, InversionOutf);
		 //InversionOutf.close();
       for (unsigned int i = 0; i < NumBoxes; i++) Inv_NT[i].clear();

} // Analyze_TD_INV_LI_Others

		 cout << "Searching short insertions ... " << endl;
		 for (unsigned ReadIndex = 0; ReadIndex < Reads.size(); ReadIndex++) {
			 if (Reads[ReadIndex].Used || Reads[ReadIndex].UP_Far.empty()) continue;
			 if (!Reads[ReadIndex].UP_Far.empty())
			 {
				 if (Reads[ReadIndex].MatchedD == Plus) {
				    for (short MAX_SNP_ERROR_index = 0; MAX_SNP_ERROR_index <= Reads[ReadIndex].MAX_SNP_ERROR; MAX_SNP_ERROR_index++) {
					    for (int CloseIndex = Reads[ReadIndex].UP_Close.size() - 1; CloseIndex >= 0; CloseIndex--) {
							//cout << "+" << CloseIndex << endl;
							 if (Reads[ReadIndex].Used) break;
							 if (Reads[ReadIndex].UP_Close[CloseIndex].Mismatches > MAX_SNP_ERROR_index) continue;
							 for (int FarIndex = Reads[ReadIndex].UP_Far.size() - 1; FarIndex >= 0; FarIndex--) {
								 //cout << "+" << FarIndex << endl;
								 if (Reads[ReadIndex].Used) break;
								 if (Reads[ReadIndex].UP_Far[FarIndex].Mismatches > MAX_SNP_ERROR_index) continue;
								 if (Reads[ReadIndex].UP_Far[FarIndex].Mismatches + Reads[ReadIndex].UP_Close[CloseIndex].Mismatches > MAX_SNP_ERROR_index)
									 continue;
								 if (Reads[ReadIndex].UP_Far[FarIndex].Direction == Minus) {
									if (Reads[ReadIndex].UP_Far[FarIndex].AbsLoc == Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc + 1 && Reads[ReadIndex].UP_Close[CloseIndex].LengthStr + Reads[ReadIndex].UP_Far[FarIndex].LengthStr < Reads[ReadIndex].ReadLength) {

										 Reads[ReadIndex].Left = Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc - Reads[ReadIndex].UP_Close[CloseIndex].LengthStr + 1;
										 Reads[ReadIndex].Right = Reads[ReadIndex].UP_Far[FarIndex].AbsLoc + Reads[ReadIndex].UP_Far[FarIndex].LengthStr - 1;
										 Reads[ReadIndex].BP = Reads[ReadIndex].UP_Close[CloseIndex].LengthStr - 1;

										 Reads[ReadIndex].IndelSize = Reads[ReadIndex].ReadLengthMinus - (Reads[ReadIndex].Right - Reads[ReadIndex].Left);
										 Reads[ReadIndex].InsertedStr = ReverseComplement(Reads[ReadIndex].UnmatchedSeq).substr(Reads[ReadIndex].BP + 1, Reads[ReadIndex].IndelSize);
										 Reads[ReadIndex].NT_str = "";
										 //                                 LeftReads[Left_Index].InsertedStr = "";
										 Reads[ReadIndex].BPLeft = Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc - SpacerBeforeAfter;
										 Reads[ReadIndex].BPRight = Reads[ReadIndex].UP_Far[FarIndex].AbsLoc - SpacerBeforeAfter;
										 Reads[ReadIndex].score = Const_I + LOG14 * Reads[ReadIndex].ReadLength - LOG14 * Reads[ReadIndex].IndelSize + Const_Log_T;
										 //LeftReads[Left_Index].OK = true;
                                        //cout << "here" << endl;
                                        unsigned RealBP_left= Reads[ReadIndex].BPLeft;
                                        unsigned RealBP_right = Reads[ReadIndex].BPRight;//, DIFF;
                                        GetRealStart4Insertion(CurrentChr, Reads[ReadIndex].InsertedStr, RealBP_left, RealBP_right);
                                        short DIFF = Reads[ReadIndex].BPLeft - RealBP_left;
                                        DIFF = !((Reads[ReadIndex].BP - 1)<DIFF)?DIFF:(Reads[ReadIndex].BP - 1); // min(DIFF, currentRead.BP - 1);
                                        //cout << DIFF << endl;
                                        if (DIFF) {
                                            //std::cout << DIFF << std::endl;
                                            Reads[ReadIndex].BP -= DIFF;
                                            Reads[ReadIndex].BPLeft -= DIFF;
                                            Reads[ReadIndex].BPRight  -= DIFF;
                                            Reads[ReadIndex].InsertedStr = ReverseComplement(Reads[ReadIndex].UnmatchedSeq).substr(Reads[ReadIndex].BP + 1, Reads[ReadIndex].IndelSize);
                                        }
                                        //cout << "there" << endl;
										 {
											 //if (RangeIndex == 0)
											 {
												 SIs[Reads[ReadIndex].BPLeft / BoxSize].push_back(ReadIndex);
												 Reads[ReadIndex].Used = true;
												 Count_SI_Plus++;
												 Count_SI++;
											 }
										 }
									 }
								}
							 }
						 }
					 }
				 }
				 else if (Reads[ReadIndex].MatchedD == Minus) {
					 for (short MAX_SNP_ERROR_index = 0; MAX_SNP_ERROR_index <= Reads[ReadIndex].MAX_SNP_ERROR; MAX_SNP_ERROR_index++) {
						 for (int CloseIndex = Reads[ReadIndex].UP_Close.size() - 1; CloseIndex >= 0; CloseIndex--) {
						    //cout << "-: " << CloseIndex << endl;
						    if (Reads[ReadIndex].Used) break;
							 if (Reads[ReadIndex].UP_Close[CloseIndex].Mismatches > MAX_SNP_ERROR_index) continue;
						    for (int FarIndex = Reads[ReadIndex].UP_Far.size() - 1; FarIndex >= 0; FarIndex--) {
							    //cout << "-: " << FarIndex << endl;
							    if (Reads[ReadIndex].Used) break;
								 if (Reads[ReadIndex].UP_Far[FarIndex].Mismatches > MAX_SNP_ERROR_index) continue;
							    if (Reads[ReadIndex].UP_Far[FarIndex].Mismatches + Reads[ReadIndex].UP_Close[CloseIndex].Mismatches > MAX_SNP_ERROR_index)
								    continue;
							    if (Reads[ReadIndex].UP_Far[FarIndex].Direction == Plus) {
								    if (Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc == Reads[ReadIndex].UP_Far[FarIndex].AbsLoc + 1 && Reads[ReadIndex].UP_Far[FarIndex].LengthStr + Reads[ReadIndex].UP_Close[CloseIndex].LengthStr < Reads[ReadIndex].ReadLength) {

									    Reads[ReadIndex].Left = Reads[ReadIndex].UP_Far[FarIndex].AbsLoc - Reads[ReadIndex].UP_Far[FarIndex].LengthStr + 1;
									    Reads[ReadIndex].Right = Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc + Reads[ReadIndex].UP_Close[CloseIndex].LengthStr - 1;
									    Reads[ReadIndex].BP = Reads[ReadIndex].UP_Far[FarIndex].LengthStr - 1;

									    Reads[ReadIndex].IndelSize = Reads[ReadIndex].ReadLengthMinus - (Reads[ReadIndex].Right - Reads[ReadIndex].Left);
									    Reads[ReadIndex].InsertedStr = Reads[ReadIndex].UnmatchedSeq.substr(Reads[ReadIndex].BP + 1, Reads[ReadIndex].IndelSize);
									    Reads[ReadIndex].NT_str = "";
									    //   LeftReads[Left_Index].InsertedStr = "";
									    Reads[ReadIndex].BPLeft = Reads[ReadIndex].UP_Far[FarIndex].AbsLoc - SpacerBeforeAfter;// - Fragments[LeftReads[Left_Index].MatchedSeqID].Start;
									    Reads[ReadIndex].BPRight = Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc - SpacerBeforeAfter;// - Fragments[LeftReads[Left_Index].MatchedSeqID].Start;
									    Reads[ReadIndex].score = Const_I + LOG14 * Reads[ReadIndex].ReadLength - LOG14 * Reads[ReadIndex].IndelSize + Const_Log_T;
									    //LeftReads[Left_Index].OK = true;
                                        //cout << "here" << endl;
                                        unsigned RealBP_left= Reads[ReadIndex].BPLeft;
                                        unsigned RealBP_right = Reads[ReadIndex].BPRight;//, DIFF;
                                        //if (Reads[ReadIndex].InsertedStr.size()) {
                                            GetRealStart4Insertion(CurrentChr, Reads[ReadIndex].InsertedStr, RealBP_left, RealBP_right);
                                        //}
                                        //else {
                                        //    GetRealStart4Deletion(CurrentChr, RealBP_left, RealBP_right);
                                        //}
                                        short DIFF = Reads[ReadIndex].BPLeft - RealBP_left;
                                        DIFF = !((Reads[ReadIndex].BP - 1)<DIFF)?DIFF:(Reads[ReadIndex].BP - 1); // min(DIFF, currentRead.BP - 1);
                                        //cout << DIFF << endl;
                                        if (DIFF) {
                                            //std::cout << DIFF << std::endl;
                                            Reads[ReadIndex].BP -= DIFF;
                                            Reads[ReadIndex].BPLeft -= DIFF;
                                            Reads[ReadIndex].BPRight  -= DIFF;
                                            Reads[ReadIndex].InsertedStr = Reads[ReadIndex].UnmatchedSeq.substr(Reads[ReadIndex].BP + 1, Reads[ReadIndex].IndelSize);
                                        }
                                        //cout << "there" << endl;
									    {
										    SIs[Reads[ReadIndex].BPLeft / BoxSize].push_back(ReadIndex);
										    Reads[ReadIndex].Used = true;
										    Count_SI++;
										    Count_SI_Minus++;
									    }
								    }
							    }
						    }
					    }
					 }
				 }
			 }
		 }
		 cout << "Total: " << Count_SI << "\t+" << Count_SI_Plus << "\t-" << Count_SI_Minus << endl;
		 ofstream SIoutputfile(SIOutputFilename.c_str());
		 SortOutputSI(NumBoxes, CurrentChr, Reads, SIs, SIoutputfile);
		 SIoutputfile.close();
		 for (unsigned int i = 0; i < NumBoxes; i++) SIs[i].clear();

		 unsigned Count_Close = 0;
		 unsigned Count_Far = 0;
		 unsigned Count_Used = 0;
		 vector <SPLIT_READ> UnusedClosedReads;
		 for (unsigned Index = 0; Index < Reads.size(); Index++) {
			 if (Reads[Index].UP_Close.size()) Count_Close++;
			 if (Reads[Index].UP_Far.size()) Count_Far++;
			 else UnusedClosedReads.push_back(Reads[Index]);
			 if (Reads[Index].Used) Count_Used++;
		 }

		 cout << "Total: " << Reads.size() << "\tClose end found " << Count_Close << "\tFar end found " << Count_Far << "\tUsed\t" << Count_Used << "\n\n";

		 Reads.clear();

		 ofstream LargeInsertionOutf(LargeInsertionOutputFilename.c_str());
		 SortOutputLI(CurrentChr, UnusedClosedReads, LargeInsertionOutf);
		 LargeInsertionOutf.close();

		 ofstream RestOutf(RestOutputFilename.c_str());
		 vector <SPLIT_READ> RestReads;
		 for (int i = 0; i < UnusedClosedReads.size(); i++) {
			 if (UnusedClosedReads[i].Used == false) RestReads.push_back(UnusedClosedReads[i]);
		 }
		 UnusedClosedReads.clear();
		 SortOutputRest(CurrentChr, RestReads, RestOutf);
		 RestOutf.close();



		 Time_Sort_E = time(NULL);
		 //Time_Load_E = time(NULL);
		 //cout << "SI: " << Count_SI_Plus << " " << Count_SI_Minus << "\n"
		 //<< "D: " << Count_D_Plus << " " << Count_D_Minus << endl;
		 //cout << Entering_D_Plus << " " << Entering_D_Minus << endl;
		 //cout << Plus_Sum_Left << " " << Plus_Sum_Right << "\n" << Minus_Sum_Right << " " << Minus_Sum_Left << endl;

      AllLoadings += (unsigned int)difftime(Time_Load_E, Time_Load_S);
      //AllMinings += (unsigned int)difftime(Time_Mine_E, Time_Load_E);
      AllSortReport += (unsigned int)difftime(Time_Sort_E, Time_Load_E);
      //Time_Load_S = time(NULL);
    //}

  cout << "Loading genome sequences and reads: " << AllLoadings << " seconds." << endl;
  //cout << "Mining indels: " << AllMinings << " seconds." << endl;
  cout << "Mining, Sorting and output results: " << AllSortReport << " seconds." << endl;
  return 0;
}//main


void ReadInOneChr(ifstream & inf_Seq, string & TheInput, const string & ChrName) {
   TheInput.clear();
   char TempChar;
   string TempLine, TempChrName;
	inf_Seq >> TempChar;
	if (TempChar != '>') {
		cout << "Please use fasta format for the reference file." << endl;
		TheInput.clear();
		return;
	}
	while (inf_Seq >> TempChrName) {
		cout << "Processing chromosome " << TempChrName << endl;
		if (!TheInput.empty()) {
		   cout << "Skip the rest of chromosomes.\n";
			break;
		}
		if (TempChrName == ChrName) {
			getline(inf_Seq, TempLine);
			while (inf_Seq >> TempChar) {
				if (TempChar != '\n' && TempChar != '\r') {
					if (TempChar == '>') {
						break;
					}
					else {
						if ( 'a' <= TempChar && TempChar <= 'z' )
							TempChar = TempChar + ( 'A' - 'a' );
						switch (TempChar) {
							case 'A': TheInput += 'A'; break;  // 00000000
							case 'C': TheInput += 'C'; break; // 00010000
							case 'G': TheInput += 'G'; break; // 00100000
							case 'T': TheInput += 'T'; break; // 00110000
							default:  TheInput += 'N';        // 01000000
						}
					} // else TempChar
				}
			}
		}
		else {
			getline(inf_Seq, TempLine);
			while (inf_Seq >> TempChar) {
				if (TempChar != '\n' && TempChar != '\r') {
					if (TempChar == '>') {
						break;
					}
					else {

					} // else TempChar
				}
			}
		}
	}
   cout << ChrName << "\t" << TheInput.size() << "\t";
	if (!TheInput.empty()) {
		string Spacer = "";
		for (unsigned i = 0; i < SpacerBeforeAfter; i++) Spacer += "N";
		TheInput = Spacer + TheInput + Spacer;
	}
	cout << TheInput.size() << endl;
   return;
}

short ReadInRead(ifstream & inf_ReadSeq, const string & FragName, const string & CurrentChr, vector <SPLIT_READ> & Reads) {
	cout << "Scanning and processing reads anchored in " << FragName << endl;
	//short ADDITIONAL_MISMATCH = 1;

	SPLIT_READ Temp_One_Read;
	unsigned int NumReadScanned = 0;
	unsigned int NumReadInChr = 0;
	unsigned int InChrPlus = 0;
	unsigned int InChrMinus = 0;
	unsigned int GetPlus = 0;
	unsigned int GetMinus = 0;
	ReportLength = 0;
	//cout << LeftReads.size() << endl;
	string TempQC, TempLine, TempStr, TempFragName;
	//int TempInt;
	VectorTag.clear();
	while (inf_ReadSeq >> Temp_One_Read.Name) {
		if (Temp_One_Read.Name[0] != FirstCharReadName) {
			cout << "Something wrong with the read name: " << Temp_One_Read.Name << endl;
			Reads.clear();
			return 1;
		}
		NumReadScanned++;
		if (NumReadScanned % 1000000 == 0)
			cout << NumReadScanned << endl;
		getline(inf_ReadSeq, TempLine);
		inf_ReadSeq >> Temp_One_Read.UnmatchedSeq;
		getline(inf_ReadSeq, TempLine);
		//inf_ReadSeq >> TempQC;
		//getline(inf_ReadSeq, TempLine);
		inf_ReadSeq >> Temp_One_Read.MatchedD;
		if (Temp_One_Read.MatchedD != Minus && Temp_One_Read.MatchedD != Plus) {
			cout << Temp_One_Read.Name << "\n"
			<< Temp_One_Read.UnmatchedSeq << "\n"
			<< Temp_One_Read.MatchedD << " ...\n";
			cout << "+/-" << endl;
			return 1;
		}
		//   >> TempInt
		inf_ReadSeq >> Temp_One_Read.FragName
		>> Temp_One_Read.MatchedRelPos
		>> Temp_One_Read.MS
		>> Temp_One_Read.InsertSize
		>> Temp_One_Read.Tag;
		//Temp_One_Read.Tag = Temp_One_Read.Tag.substr(0, 7);
		//Temp_One_Read.MS = 1;
		//Temp_One_Read.InsertSize = 300;
		getline(inf_ReadSeq, TempLine);


		if (Temp_One_Read.FragName == FragName/* && Temp_One_Read.MatchedRelPos > 10000000 && Temp_One_Read.MatchedRelPos < 15000000*/) {
			NumReadInChr++;
			Temp_One_Read.MAX_SNP_ERROR = (short)(Temp_One_Read.UnmatchedSeq.size() * Seq_Error_Rate);
			//MAX_SNP_ERROR = (short)(Temp_One_Read.UnmatchedSeq.size() * Seq_Error_Rate);

			Temp_One_Read.TOTAL_SNP_ERROR_CHECKED = Temp_One_Read.MAX_SNP_ERROR + ADDITIONAL_MISMATCH + 1;
			Temp_One_Read.TOTAL_SNP_ERROR_CHECKED_Minus = Temp_One_Read.MAX_SNP_ERROR + ADDITIONAL_MISMATCH;
			Temp_One_Read.MinClose = short(log((double)(Temp_One_Read.InsertSize * 3))/log(4.0) + 0.8) + 3;// + MAX_SNP_ERROR;//atoi(argv[1]);
			//MinFar_I = MinClose + 1;//atoi(argv[2]);
			if (Temp_One_Read.MatchedD == Plus)  {
				InChrPlus++;
				//Temp_One_Read.UnmatchedSeq = ReverseComplement(Temp_One_Read.UnmatchedSeq);
			}
			else InChrMinus++;
			if (Temp_One_Read.MatchedRelPos > CONS_Chr_Size) Temp_One_Read.MatchedRelPos = CONS_Chr_Size;
			if (Temp_One_Read.MatchedRelPos < 1) Temp_One_Read.MatchedRelPos = 0;
			GetCloseEnd(CurrentChr, Temp_One_Read);
			if (Temp_One_Read.UP_Close.size()) {
				if (ReportLength < Temp_One_Read.ReadLength) ReportLength = Temp_One_Read.ReadLength;
				Temp_One_Read.Used = false;
				Temp_One_Read.OK = true;
				//cout << Temp_One_Read.MatchedD << "\t" << Temp_One_Read.UP_Close.size() << "\t";
				CleanUniquePoints(Temp_One_Read.UP_Close);
				//cout << Temp_One_Read.UP_Close.size() << "\t" << Temp_One_Read.UP_Close[0].Direction << endl;
				Temp_One_Read.CloseEndLength =
				Temp_One_Read.UP_Close[Temp_One_Read.UP_Close.size() - 1].LengthStr;

				//if (Temp_One_Read.UP_Close.size()) {
					Reads.push_back(Temp_One_Read);
					if (Temp_One_Read.MatchedD == Plus) GetPlus++;
					else GetMinus++;
				//}
				//if (NumReadScanned % 1000 == 0)
				//cout << "NumReadScanned: " << NumReadScanned << "\n"
				//     << "NumReadInChr: " << NumReadInChr << "\n"
				//     << "NumReadStored: " << Reads.size() << "\n"
				//     << InChrPlus << "\t" << InChrMinus << "\t"
				//	  << GetPlus << "\t" << GetMinus << "\t" << GetPlus * 2.0 / NumReadInChr
				//     << "\t" << GetMinus * 2.0 / NumReadInChr
				//     << "\t" << Reads.size() * 1.0 / NumReadInChr << endl;
			}
			if (NotInVector(Temp_One_Read.Tag, VectorTag)) {
				VectorTag.push_back(Temp_One_Read.Tag);
			}
		}
	}
	if (FirstChr) {
		cout << "\nThe last read Pindel scanned: \n";
		cout << Temp_One_Read.Name << "\n"
		<< Temp_One_Read.UnmatchedSeq << "\n"
		<< Temp_One_Read.MatchedD << "\t"
		<< Temp_One_Read.FragName << "\t"
		<< Temp_One_Read.MatchedRelPos << "\t"
		<< Temp_One_Read.MS << "\t"
		<< Temp_One_Read.InsertSize << "\t"
		<< Temp_One_Read.Tag << "\n\n";
		FirstChr = false;
      //ReportLength = Temp_One_Read.UnmatchedSeq.size();
	}
	cout << "NumReadScanned:\t" << NumReadScanned << endl;
	cout << "NumReadInChr:\t" << NumReadInChr << endl;
	cout << "NumReadStored:\t" << Reads.size() << endl;
	cout << "NumReadStored / NumReadInChr = " << Reads.size() * 100.0 / NumReadInChr << " %\n"
	<< "InChrPlus \t" << InChrPlus << "\tGetPlus \t" << GetPlus << "\t" << GetPlus * 100.0 / InChrPlus
	<< " %\n" << "InChrMinus\t" << InChrMinus  << "\tGetMinus\t" << GetMinus
	<< "\t" << GetMinus * 100.0 / InChrMinus << " %\n" << endl;
	inf_ReadSeq.close();
	if (Reads.size() == 0) return 0;
	//cout << LeftReads.size() << endl;
	cout << "sorting tags ... ";
	string Str4Exchange;
	for (short i = 0; i < VectorTag.size() - 1; i++) {
		for (short j = 1; j < VectorTag.size(); j++) {
			if (VectorTag[i].size() > VectorTag[j].size()) {
				Str4Exchange = VectorTag[i];
				VectorTag[i] = VectorTag[j];
				VectorTag[j] = Str4Exchange;
			}
			else if (VectorTag[i].size() == VectorTag[j].size()) {
				for (short k = 0; k < VectorTag[i].size(); k++) {
					if ((short)VectorTag[i][k] > (short)VectorTag[j][k]) {
       	         Str4Exchange = VectorTag[i];
						VectorTag[i] = VectorTag[j];
						VectorTag[j] = Str4Exchange;
						break;
					}
				}
			}
		}
	}
	cout << " finished!" << endl;
	return 0;
}



vector <string> ReverseComplement(const vector <string> & InputPatterns) {
   vector <string> OutputPatterns;// = InputPatterns;
   unsigned int NumPattern = InputPatterns.size();
   OutputPatterns.resize(NumPattern);
   for (unsigned int i = 0; i < NumPattern; i++) {
      OutputPatterns[i] = ReverseComplement(InputPatterns[i]);
   }
   return OutputPatterns;
}


string Reverse(const string & InputPattern) {
   string OutputPattern = InputPattern;
   unsigned int LenPattern = InputPattern.size();
   for (unsigned int j = 0; j < LenPattern; j++)
      OutputPattern[j] = InputPattern[j];
   return OutputPattern;
}

string ReverseComplement(const string & InputPattern) {
   string OutputPattern = InputPattern;
   unsigned int LenPattern = InputPattern.size();
   for (unsigned int j = 0; j < LenPattern; j++)
      OutputPattern[j] = Convert2RC4N[(unsigned int)InputPattern[LenPattern - j - 1]];
   return OutputPattern;
}

void CheckLeft_Close(const SPLIT_READ & OneRead,
							const string & TheInput,
							const string & CurrentReadSeq,
							const vector <unsigned int> Left_PD[],
							const short & BP_Left_Start,
							const short & BP_Left_End,
							const short & CurrentLength,
							vector <UniquePoint> & LeftUP) {
	int Sum;
   if (CurrentLength >= BP_Left_Start && CurrentLength <= BP_Left_End) {
      // put it to LeftUP if unique
		for (short i = 0; i <= OneRead.MAX_SNP_ERROR; i++) {
			if (Left_PD[i].size() == 1 && CurrentLength >= BP_Left_Start + i) {
				Sum = 0;
				if (ADDITIONAL_MISMATCH)
		   		for (short j = 1; j <= ADDITIONAL_MISMATCH; j++)
			   		Sum += Left_PD[i + j].size();

				if (Sum == 0 && i <= (short)(CurrentLength * Seq_Error_Rate + 1)) {
					UniquePoint TempOne;
					TempOne.LengthStr = CurrentLength;
					TempOne.AbsLoc = Left_PD[i][0];
					TempOne.Direction = FORWARD;
					TempOne.Strand = ANTISENSE;
					TempOne.Mismatches = i;
					LeftUP.push_back(TempOne);
					break;
				}
			}
		}
   }
   if (CurrentLength < BP_Left_End) {
      vector <unsigned int> Left_PD_Output[OneRead.TOTAL_SNP_ERROR_CHECKED];
		for (int CheckedIndex = 0; CheckedIndex < OneRead.TOTAL_SNP_ERROR_CHECKED; CheckedIndex++) {
			Left_PD_Output[CheckedIndex].reserve(Left_PD[CheckedIndex].size());
		}
      const char CurrentChar = CurrentReadSeq[CurrentLength];
      //const int SizeOfCurrent = Left_PD.size();
      unsigned int pos;
		for (int i = 0; i < OneRead.TOTAL_SNP_ERROR_CHECKED_Minus; i++) {
			int SizeOfCurrent = Left_PD[i].size();
			if (CurrentChar == 'N') {
				for (int j = 0; j < SizeOfCurrent; j++) {
					pos =  Left_PD[i][j] + 1;
					if (Match2N[(short)TheInput[pos]] == 'N')
						Left_PD_Output[i].push_back(pos);
				}
			}
			else {
				for (int j = 0; j < SizeOfCurrent; j++) {
					pos =  Left_PD[i][j] + 1;
					if (TheInput[pos] == CurrentChar)
						Left_PD_Output[i].push_back(pos);
					else Left_PD_Output[i + 1].push_back(pos);
				}
			}
		}

		int SizeOfCurrent = Left_PD[OneRead.TOTAL_SNP_ERROR_CHECKED_Minus].size();
		if (CurrentChar == 'N') {
			for (int j = 0; j < SizeOfCurrent; j++) {
				pos =  Left_PD[OneRead.TOTAL_SNP_ERROR_CHECKED_Minus][j] + 1;
				if (Match2N[(short)TheInput[pos]] == 'N')
					Left_PD_Output[OneRead.TOTAL_SNP_ERROR_CHECKED_Minus].push_back(pos);
			}
		}
		else {
			for (int j = 0; j < SizeOfCurrent; j++) {
				pos =  Left_PD[OneRead.TOTAL_SNP_ERROR_CHECKED_Minus][j] + 1;
				if (TheInput[pos] == CurrentChar)
					Left_PD_Output[OneRead.TOTAL_SNP_ERROR_CHECKED_Minus].push_back(pos);
				//else Left_PD_Output[i + 1].push_back(pos);
			}
		}
		Sum = 0;
		for (int i = 0; i <= OneRead.MAX_SNP_ERROR; i++) {
			Sum += Left_PD_Output[i].size();
		}
      if (Sum) {
         const short CurrentLengthOutput = CurrentLength + 1;
         CheckLeft_Close(OneRead, TheInput, CurrentReadSeq, Left_PD_Output,
								 BP_Left_Start, BP_Left_End,
								 CurrentLengthOutput, LeftUP);
      }
      else return;
   }
   else return;
}

void CheckRight_Close(const SPLIT_READ & OneRead,
							 const string & TheInput,
							 const string & CurrentReadSeq,
							 const vector <unsigned int> Right_PD[],
							 const short & BP_Right_Start,
							 const short & BP_Right_End,
							 const short & CurrentLength,
							 vector <UniquePoint> & RightUP) {
	//cout << CurrentLength << "\t" << RightUP.size() << "\t" << Right_PD[0].size() << "\t" << Right_PD[1].size() << endl;
   short ReadLengthMinus = CurrentReadSeq.size() - 1;
	int Sum;
   if (CurrentLength >= BP_Right_Start && CurrentLength <= BP_Right_End) {
	   for (short i = 0; i <= OneRead.MAX_SNP_ERROR; i++) {
		   if (Right_PD[i].size() == 1 && CurrentLength >= BP_Right_Start + i) {
				Sum = 0;
				if (ADDITIONAL_MISMATCH)
					for (short j = 1; j <= ADDITIONAL_MISMATCH; j++)
						Sum += Right_PD[i + j].size();

				if (Sum == 0 && i <= (short)(CurrentLength * Seq_Error_Rate + 1)) {
					UniquePoint TempOne;
					TempOne.LengthStr = CurrentLength;
					TempOne.AbsLoc = Right_PD[i][0];
					TempOne.Direction = BACKWARD;
					TempOne.Strand = SENSE;
					TempOne.Mismatches = i;
					RightUP.push_back(TempOne);
					break;
				}
			}
		}
   }

   if (CurrentLength < BP_Right_End) {
      vector <unsigned int> Right_PD_Output[OneRead.TOTAL_SNP_ERROR_CHECKED];
		for (int CheckedIndex = 0; CheckedIndex < OneRead.TOTAL_SNP_ERROR_CHECKED; CheckedIndex++) {
			Right_PD_Output[CheckedIndex].reserve(Right_PD[CheckedIndex].size());
		}
      const char CurrentChar = CurrentReadSeq[ReadLengthMinus - CurrentLength];
      unsigned int pos;

		for (int i = 0; i < OneRead.TOTAL_SNP_ERROR_CHECKED_Minus; i++) {
			int SizeOfCurrent = Right_PD[i].size();
			if (CurrentChar == 'N') {
				for (int j = 0; j < SizeOfCurrent; j++) {
					pos =  Right_PD[i][j] - 1;
					if (Match2N[(short)TheInput[pos]] == 'N')
						Right_PD_Output[i].push_back(pos);
				}
			}
			else {
				for (int j = 0; j < SizeOfCurrent; j++) {
					pos =  Right_PD[i][j] - 1;
					if (TheInput[pos] == CurrentChar)
						Right_PD_Output[i].push_back(pos);
					else Right_PD_Output[i + 1].push_back(pos);
				}
			}
		}

		int SizeOfCurrent = Right_PD[OneRead.TOTAL_SNP_ERROR_CHECKED_Minus].size();
		if (CurrentChar == 'N') {
			for (int j = 0; j < SizeOfCurrent; j++) {
				pos =  Right_PD[OneRead.TOTAL_SNP_ERROR_CHECKED_Minus][j] - 1;
				if (Match2N[(short)TheInput[pos]] == 'N')
					Right_PD_Output[OneRead.TOTAL_SNP_ERROR_CHECKED_Minus].push_back(pos);
			}
		}
		else {
			for (int j = 0; j < SizeOfCurrent; j++) {
				pos =  Right_PD[OneRead.TOTAL_SNP_ERROR_CHECKED_Minus][j] - 1;
				if (TheInput[pos] == CurrentChar)
					Right_PD_Output[OneRead.TOTAL_SNP_ERROR_CHECKED_Minus].push_back(pos);
				//else Left_PD_Output[i + 1].push_back(pos);
			}
		}

		Sum = 0;
		for (int i = 0; i <= OneRead.MAX_SNP_ERROR; i++) {
			Sum += Right_PD_Output[i].size();
		}
      if (Sum) {
         short CurrentLength_output = CurrentLength + 1;
         CheckRight_Close(OneRead, TheInput, CurrentReadSeq, Right_PD_Output,
								  BP_Right_Start, BP_Right_End,
								  CurrentLength_output, RightUP);
      }
      else return;
   }
   else return;
}

void CheckLeft_Far(const SPLIT_READ & OneRead,
						 const string & TheInput,
							const string & CurrentReadSeq,
							const vector <unsigned int> Left_PD[],
							const short & BP_Left_Start,
							const short & BP_Left_End,
							const short & CurrentLength,
							vector <UniquePoint> & LeftUP) {
	int Sum;
   if (CurrentLength >= BP_Left_Start && CurrentLength <= BP_Left_End) {
      // put it to LeftUP if unique
		for (short i = 0; i <= OneRead.MAX_SNP_ERROR; i++) {
			if (Left_PD[i].size() == 1 && CurrentLength >= BP_Left_Start + i) {
				Sum = 0;
				if (ADDITIONAL_MISMATCH)
		   		for (short j = 1; j <= ADDITIONAL_MISMATCH; j++)
			   		Sum += Left_PD[i + j].size();

				if (Sum == 0 && i <= (short)(CurrentLength * Seq_Error_Rate + 1)) {
					UniquePoint TempOne;
					TempOne.LengthStr = CurrentLength;
					TempOne.AbsLoc = Left_PD[i][0];
					TempOne.Direction = FORWARD;
					TempOne.Strand = SENSE;
					TempOne.Mismatches = i;
					LeftUP.push_back(TempOne);
					break;
				}
			}
		}
   }
   if (CurrentLength < BP_Left_End) {
      vector <unsigned int> Left_PD_Output[OneRead.TOTAL_SNP_ERROR_CHECKED];
		for (int CheckedIndex = 0; CheckedIndex < OneRead.TOTAL_SNP_ERROR_CHECKED; CheckedIndex++) {
			Left_PD_Output[CheckedIndex].reserve(Left_PD[CheckedIndex].size());
		}
      const char CurrentChar = CurrentReadSeq[CurrentLength];
      //const int SizeOfCurrent = Left_PD.size();
      unsigned int pos;
		//if (TOTAL_SNP_ERROR_CHECKED_Minus)
		{
			for (int i = 0; i < OneRead.TOTAL_SNP_ERROR_CHECKED_Minus; i++) {
				int SizeOfCurrent = Left_PD[i].size();
				if (CurrentChar == 'N') {
					for (int j = 0; j < SizeOfCurrent; j++) {
						pos =  Left_PD[i][j] + 1;
						if (Match2N[(short)TheInput[pos]] == 'N')
							Left_PD_Output[i].push_back(pos);
					}
				}
				else {
					for (int j = 0; j < SizeOfCurrent; j++) {
						pos =  Left_PD[i][j] + 1;
						if (TheInput[pos] == CurrentChar)
							Left_PD_Output[i].push_back(pos);
						else Left_PD_Output[i + 1].push_back(pos);
					}
				}
			}

			int SizeOfCurrent = Left_PD[OneRead.TOTAL_SNP_ERROR_CHECKED_Minus].size();
			if (CurrentChar == 'N') {
				for (int j = 0; j < SizeOfCurrent; j++) {
					pos =  Left_PD[OneRead.TOTAL_SNP_ERROR_CHECKED_Minus][j] + 1;
					if (Match2N[(short)TheInput[pos]] == 'N')
						Left_PD_Output[OneRead.TOTAL_SNP_ERROR_CHECKED_Minus].push_back(pos);
				}
			}
			else {
				for (int j = 0; j < SizeOfCurrent; j++) {
					pos =  Left_PD[OneRead.TOTAL_SNP_ERROR_CHECKED_Minus][j] + 1;
					if (TheInput[pos] == CurrentChar)
						Left_PD_Output[OneRead.TOTAL_SNP_ERROR_CHECKED_Minus].push_back(pos);
					//else Left_PD_Output[i + 1].push_back(pos);
				}
			}

			Sum = 0;
			for (int i = 0; i <= OneRead.MAX_SNP_ERROR; i++) {
				Sum += Left_PD_Output[i].size();
			}
			if (Sum) {
				const short CurrentLengthOutput = CurrentLength + 1;
				CheckLeft_Far(OneRead,TheInput, CurrentReadSeq, Left_PD_Output,
								  BP_Left_Start, BP_Left_End,
								  CurrentLengthOutput, LeftUP);
			}
			else return;
		}
		/*
      else { // TOTAL_SNP_ERROR_CHECKED_Minus

			int SizeOfCurrent = Left_PD[TOTAL_SNP_ERROR_CHECKED_Minus].size();
			if (CurrentChar == 'N') {
				for (int j = 0; j < SizeOfCurrent; j++) {
					pos =  Left_PD[TOTAL_SNP_ERROR_CHECKED_Minus][j] + 1;
					if (Match2N[(short)TheInput[pos]] == 'N')
						Left_PD_Output[TOTAL_SNP_ERROR_CHECKED_Minus].push_back(pos);
				}
			}
			else {
				for (int j = 0; j < SizeOfCurrent; j++) {
					pos =  Left_PD[TOTAL_SNP_ERROR_CHECKED_Minus][j] + 1;
					if (TheInput[pos] == CurrentChar)
						Left_PD_Output[TOTAL_SNP_ERROR_CHECKED_Minus].push_back(pos);
					//else Left_PD_Output[i + 1].push_back(pos);
				}
			}
			if (!Left_PD_Output[0].empty()) {
				const short CurrentLengthOutput = CurrentLength + 1;
				CheckLeft_Far(TheInput, CurrentReadSeq, Left_PD_Output,
								  BP_Left_Start, BP_Left_End,
								  CurrentLengthOutput, LeftUP);
			}
			else return;
		}
		 */
   }
   else return;
}

void CheckRight_Far(const SPLIT_READ & OneRead,
						  const string & TheInput,
							 const string & CurrentReadSeq,
							 const vector <unsigned int> Right_PD[],
							 const short & BP_Right_Start,
							 const short & BP_Right_End,
							 const short & CurrentLength,
							 vector <UniquePoint> & RightUP) {
   short ReadLengthMinus = CurrentReadSeq.size() - 1;
	int Sum;
   if (CurrentLength >= BP_Right_Start && CurrentLength <= BP_Right_End) {
	   for (short i = 0; i <= OneRead.MAX_SNP_ERROR; i++) {
		   if (Right_PD[i].size() == 1 && CurrentLength >= BP_Right_Start + i) {
				Sum = 0;
				if (ADDITIONAL_MISMATCH)
					for (short j = 1; j <= ADDITIONAL_MISMATCH; j++)
						Sum += Right_PD[i + j].size();

				if (Sum == 0 && i <= (short)(CurrentLength * Seq_Error_Rate + 1)) {
					UniquePoint TempOne;
					TempOne.LengthStr = CurrentLength;
					TempOne.AbsLoc = Right_PD[i][0];
					TempOne.Direction = BACKWARD;
					TempOne.Strand = ANTISENSE;
					TempOne.Mismatches = i;
					RightUP.push_back(TempOne);
					break;
				}
			}
		}
   }

   if (CurrentLength < BP_Right_End) {
      vector <unsigned int> Right_PD_Output[OneRead.TOTAL_SNP_ERROR_CHECKED];
		for (int CheckedIndex = 0; CheckedIndex < OneRead.TOTAL_SNP_ERROR_CHECKED; CheckedIndex++) {
			Right_PD_Output[CheckedIndex].reserve(Right_PD[CheckedIndex].size());
		}
      const char CurrentChar = CurrentReadSeq[ReadLengthMinus - CurrentLength];
      unsigned int pos;

		//if (TOTAL_SNP_ERROR_CHECKED_Minus)
		{
			for (int i = 0; i < OneRead.TOTAL_SNP_ERROR_CHECKED_Minus; i++) {
				int SizeOfCurrent = Right_PD[i].size();
				if (CurrentChar == 'N') {
					for (int j = 0; j < SizeOfCurrent; j++) {
						pos =  Right_PD[i][j] - 1;
						if (Match2N[(short)TheInput[pos]] == 'N')
							Right_PD_Output[i].push_back(pos);
					}
				}
				else {
					for (int j = 0; j < SizeOfCurrent; j++) {
						pos =  Right_PD[i][j] - 1;
						if (TheInput[pos] == CurrentChar)
							Right_PD_Output[i].push_back(pos);
						else Right_PD_Output[i + 1].push_back(pos);
					}
				}
			}

			int SizeOfCurrent = Right_PD[OneRead.TOTAL_SNP_ERROR_CHECKED_Minus].size();
			if (CurrentChar == 'N') {
				for (int j = 0; j < SizeOfCurrent; j++) {
					pos =  Right_PD[OneRead.TOTAL_SNP_ERROR_CHECKED_Minus][j] - 1;
					if (Match2N[(short)TheInput[pos]] == 'N')
						Right_PD_Output[OneRead.TOTAL_SNP_ERROR_CHECKED_Minus].push_back(pos);
				}
			}
			else {
				for (int j = 0; j < SizeOfCurrent; j++) {
					pos =  Right_PD[OneRead.TOTAL_SNP_ERROR_CHECKED_Minus][j] - 1;
					if (TheInput[pos] == CurrentChar)
						Right_PD_Output[OneRead.TOTAL_SNP_ERROR_CHECKED_Minus].push_back(pos);
					//else Left_PD_Output[i + 1].push_back(pos);
				}
			}

			Sum = 0;
			for (int i = 0; i <= OneRead.MAX_SNP_ERROR; i++) {
				Sum += Right_PD_Output[i].size();
			}
			if (Sum) {
				short CurrentLength_output = CurrentLength + 1;
				CheckRight_Far(OneRead, TheInput, CurrentReadSeq, Right_PD_Output,
									BP_Right_Start, BP_Right_End,
									CurrentLength_output, RightUP);
			}
			else return;
		}
		/*
      else { // TOTAL_SNP_ERROR_CHECKED_Minus

			int SizeOfCurrent = Right_PD[TOTAL_SNP_ERROR_CHECKED_Minus].size();
			if (CurrentChar == 'N') {
				for (int j = 0; j < SizeOfCurrent; j++) {
					pos =  Right_PD[TOTAL_SNP_ERROR_CHECKED_Minus][j] - 1;
					if (Match2N[(short)TheInput[pos]] == 'N')
						Right_PD_Output[TOTAL_SNP_ERROR_CHECKED_Minus].push_back(pos);
				}
			}
			else {
				for (int j = 0; j < SizeOfCurrent; j++) {
					pos =  Right_PD[TOTAL_SNP_ERROR_CHECKED_Minus][j] - 1;
					if (TheInput[pos] == CurrentChar)
						Right_PD_Output[TOTAL_SNP_ERROR_CHECKED_Minus].push_back(pos);
					//else Left_PD_Output[i + 1].push_back(pos);
				}
			}

			if (!Right_PD_Output[TOTAL_SNP_ERROR_CHECKED_Minus].empty()) {
				short CurrentLength_output = CurrentLength + 1;
				CheckRight_Far(TheInput, CurrentReadSeq, Right_PD_Output,
									BP_Right_Start, BP_Right_End,
									CurrentLength_output, RightUP);
			}
			else return;
		}
		 */
   }
   else return;
}


short CompareTwoReads(const SPLIT_READ & First, const SPLIT_READ & Second) {
   //if (First.MatchedSeqID < Second.MatchedSeqID) return 0;
   //else if (First.MatchedSeqID > Second.MatchedSeqID) return 1;
   short FirstDISize = First.NT_str.size();
   short SecondDISize = Second.NT_str.size();
   short FirstSISize = First.InsertedStr.size();
   short SecondSISize = Second.InsertedStr.size();

   if (First.BPLeft > Second.BPLeft) return 1;
   else if (First.BPLeft < Second.BPLeft) return 0;
   else {
      {
         if (First.IndelSize > Second.IndelSize) return 1;
         else if (First.IndelSize < Second.IndelSize) return 0;
         else
         {
            if (First.Tag.size() < Second.Tag.size()) return 0;
            else if (First.Tag.size() > Second.Tag.size()) return 1;
            else if (First.Tag.size() == Second.Tag.size()) {
               for (unsigned posindex = 0; posindex < First.Tag.size(); posindex++) {
                  if ((short)First.Tag[posindex] < (short)Second.Tag[posindex]) return 0;
                  else if ((short)First.Tag[posindex] > (short)Second.Tag[posindex]) return 1;
               }
               if (First.MatchedRelPos < Second.MatchedRelPos)  return 0;
               else if (First.MatchedRelPos > Second.MatchedRelPos) return 1;
               else {
                  if (FirstDISize > SecondDISize) return 1;
   		  else if (FirstDISize < SecondDISize) return 0;
   		  else if (FirstDISize != 0) {
                     for (int i = 0; i < FirstDISize; i++)
                        if ((short)First.NT_str[i] > (short)Second.NT_str[i]) return 1;
                        else if ((short)First.NT_str[i] < (short)Second.NT_str[i]) return 0;
                  }
                  else {
                     if (FirstSISize > SecondSISize) return 1;
                     else if (FirstSISize < SecondSISize) return 0;
                     else if (FirstSISize != 0) {
                        for (int i = 0; i < FirstSISize; i++)
                           if ((short)First.InsertedStr[i] > (short)Second.InsertedStr[i]) return 1;
                           else if ((short)First.InsertedStr[i] < (short)Second.InsertedStr[i]) return 0;
                     }
                     else return 2;
                  }
               }
               //return 2;
            }
         }
      }
   }
	return 0;
}

string Cap2Low(const string & input) {
   string output = input;
   for (unsigned int i = 0; i < output.size(); i++) {
      output[i] = Cap2LowArray[(short)input[i]];
   }
   return output;
}

bool CheckMismatches(const string & TheInput,
                     const string & InputReadSeq,
                     //const unsigned int & Start,
							const UniquePoint & UP
							) {
	//return true; 	short LengthStr;
	//unsigned int AbsLoc;
	string CurrentReadSeq;
	if (UP.Strand == SENSE) CurrentReadSeq = InputReadSeq;
	else CurrentReadSeq = ReverseComplement(InputReadSeq);
	short CurrentReadLength = CurrentReadSeq.size();
	unsigned int Start = 0;

	string BP_On_Read, BP_On_Ref;
	if (UP.Direction == FORWARD) {
		Start = UP.AbsLoc - UP.LengthStr + 1;
		BP_On_Read = CurrentReadSeq.substr(UP.LengthStr - Min_Perfect_Match_Around_BP, Min_Perfect_Match_Around_BP);
		BP_On_Ref = TheInput.substr(UP.AbsLoc - Min_Perfect_Match_Around_BP + 1, Min_Perfect_Match_Around_BP);
		if (BP_On_Read != BP_On_Ref) return false;
	}
	else if (UP.Direction == BACKWARD) {
		Start = UP.AbsLoc + UP.LengthStr - CurrentReadLength;
		BP_On_Read = CurrentReadSeq.substr(CurrentReadLength - UP.LengthStr, Min_Perfect_Match_Around_BP);
		BP_On_Ref = TheInput.substr(UP.AbsLoc, Min_Perfect_Match_Around_BP);
		if (BP_On_Read != BP_On_Ref) return false;
	}

	short MAX_ALLOWED_MISMATCHES = (short)(CurrentReadSeq.size() * MaximumAllowedMismatchRate + 1); //

   short NumMismatches = 0;                 // Match2N[(short)'A'] = 'N';

   for (short i = 0; i < CurrentReadLength; i++) {
      if (CurrentReadSeq[i] == N_char) {
         if (Match2N[(short)TheInput[Start + i]] != N_char)
            NumMismatches++;
      }
      else {
         if (TheInput[Start + i] != CurrentReadSeq[i])
            NumMismatches++;
      }
   }
   if (NumMismatches > MAX_ALLOWED_MISMATCHES) return true;
   else return false;
}

void OutputTDs(const vector <SPLIT_READ> & TDs,
                     const string & TheInput,
                     const unsigned int & C_S,
                     const unsigned int & C_E,
							const unsigned int & RealStart,
							const unsigned int & RealEnd,
                     ofstream & TDOutf) {
   //short ReadLength = Deletions[C_S].ReadLength;
   //short ReadLengthMinus = ReadLength - 1;
   unsigned int NumberOfReads = C_E - C_S + 1;
   float LeftScore = 0;
   float RightScore = 0;
   unsigned int LeftS = 1;
   unsigned int RightS = 1;
   unsigned int NumSupportPerTag[VectorTag.size()];
   for (short i = 0; i < VectorTag.size(); i++) NumSupportPerTag[i] = 0;
   for (unsigned int i = C_S; i <= C_E; i++) {
      for (short j = 0; j < VectorTag.size(); j++) {
			if (TDs[i].Tag == VectorTag[j]) NumSupportPerTag[j]++;
      }
      if (TDs[i].MatchedD == Plus) {
         LeftScore += TDs[i].score;
         LeftS++;
      }
      if (TDs[i].MatchedD == Minus) {
         RightScore += TDs[i].score;
         RightS++;
      }
   }
   unsigned int EasyScore = LeftS * RightS;
   double PreciseScore = (LeftScore + RightScore) * (-1);
   short GapSize =0;
   if (TDs[C_S].IndelSize < 14) GapSize = TDs[C_S].IndelSize;
   else GapSize = 13 + (int)log10(TDs[C_S].IndelSize - 10);

   TDOutf << "####################################################################################################" << endl;
   TDOutf << NumberOfTDInstances << "\tTD " << TDs[C_S].IndelSize// << " bases "
	<< "\tNT " << TDs[C_S].NT_size << " \"" << TDs[C_S].NT_str << "\""
	<< "\tChrID " << TDs[C_S].FragName
	<< "\tBP " << TDs[C_S].BPLeft << "\t" << TDs[C_S].BPRight + 2
	<< "\tBP_range " << TDs[C_S].BPLeft << "\t" << TDs[C_S].BPRight + 2
	<< "\tSupports " << NumberOfReads
	<< "\t+ " << LeftS - 1
	<< "\t- " << RightS - 1
	<< "\tS1 " << EasyScore << "\tS2 " << PreciseScore;

   int SUM_MS = 0;
   for (unsigned int i = C_S; i <= C_E; i++) SUM_MS += TDs[i].MS;
   TDOutf << "\tSUM_MS " << SUM_MS;

   short NumberSupSamples = 0;
   for (short i = 0; i < VectorTag.size(); ++i)
		if (NumSupportPerTag[i]) ++NumberSupSamples;
   TDOutf << "\tNumSupSamples " << NumberSupSamples;
   for (short i = 0; i < VectorTag.size(); i++)
		if (NumSupportPerTag[i])
			TDOutf << "\t" << VectorTag[i] << " " << NumSupportPerTag[i];
   TDOutf << endl;

   //DeletionOutf << TheInput.substr(Deletions[C_S].Left - ReportLength + Deletions[C_S].BP + 1, ReportLength * 2) << endl;// << endl;// ReportLength
   TDOutf << TheInput.substr(TDs[C_S].BPRight + SpacerBeforeAfter - ReportLength + 1, ReportLength);// << endl;//
	if (TDs[C_S].NT_size) {
		for (short i = 0; i < TDs[C_S].NT_size; i++) TDOutf << " ";
	}
	//ReportLength
   TDOutf << Cap2Low(TheInput.substr(TDs[C_S].BPLeft + SpacerBeforeAfter, ReportLength)) << endl;

   short SpaceBeforeReadSeq;
   for (unsigned int GoodIndex = C_S; GoodIndex <= C_E; GoodIndex++) {
      SpaceBeforeReadSeq = ReportLength - TDs[GoodIndex].BP - 1;

      for (int i = 0; i < SpaceBeforeReadSeq; i++) TDOutf << " ";
      short SpaceBeforeD = ReportLength + ReportLength - SpacerBeforeAfter - TDs[GoodIndex].ReadLength;
      if (TDs[GoodIndex].MatchedD == Minus) {
         TDOutf << TDs[GoodIndex].UnmatchedSeq << endl;
         //for (int i = 0; i < GapSize; i++) TDOutf << " ";
         //TDOutf << TDs[GoodIndex].UnmatchedSeq.substr(TDs[GoodIndex].BP + 1, TDs[GoodIndex].ReadLength - TDs[GoodIndex].BP);// << endl;
      }
      else {
         TDOutf << ReverseComplement(TDs[GoodIndex].UnmatchedSeq) << endl;
         //for (int i = 0; i < GapSize; i++) TDOutf << " ";
         //TDOutf << ReverseComplement(TDs[GoodIndex].UnmatchedSeq).substr(TDs[GoodIndex].BP + 1, TDs[GoodIndex].ReadLength - TDs[GoodIndex].BP);// << endl;
      }
      //for (int i = 0; i < SpaceBeforeD; i++) TDOutf << " ";
      TDOutf << "\t" << TDs[GoodIndex].MatchedD << "\t"
		<< TDs[GoodIndex].MatchedRelPos
		<< "\t" << TDs[GoodIndex].MS
		<< "\t" << TDs[GoodIndex].Tag
		<< "\t" <<  TDs[GoodIndex].Name << endl;
		//<< "\t" << TDs[C_S].BPLeft
		//<< "\t" << TDs[C_S].BPRight << endl;
   }
}

void OutputDeletions(const vector <SPLIT_READ> & Deletions,
                     const string & TheInput,
                     const unsigned int & C_S,
                     const unsigned int & C_E,
		     const unsigned int & RealStart,
		     const unsigned int & RealEnd,
                     ofstream & DeletionOutf) {
   //short ReadLength = Deletions[C_S].ReadLength;
   //short ReadLengthMinus = ReadLength - 1;
   unsigned int NumberOfReads = C_E - C_S + 1;
   float LeftScore = 0;
   float RightScore = 0;
   unsigned int LeftS = 1;
   unsigned int RightS = 1;
   unsigned int NumSupportPerTag[VectorTag.size()];
   for (short i = 0; i < VectorTag.size(); i++) NumSupportPerTag[i] = 0;
   for (unsigned int i = C_S; i <= C_E; i++) {
      for (short j = 0; j < VectorTag.size(); j++) {
	if (Deletions[i].Tag == VectorTag[j]) NumSupportPerTag[j]++;
      }
      if (Deletions[i].MatchedD == Plus) {
         LeftScore += Deletions[i].score;
         LeftS++;
      }
      if (Deletions[i].MatchedD == Minus) {
         RightScore += Deletions[i].score;
         RightS++;
      }
   }
   unsigned int EasyScore = LeftS * RightS;
   double PreciseScore = (LeftScore + RightScore) * (-1);
   short GapSize =0;
   if (Deletions[C_S].IndelSize < 14) GapSize = Deletions[C_S].IndelSize;
   else GapSize = 13 + (int)log10(Deletions[C_S].IndelSize - 10);

   DeletionOutf << "####################################################################################################" << endl;
   DeletionOutf << NumberOfDeletionsInstances << "\tD " << Deletions[C_S].IndelSize// << " bases "
	             << "\tNT " << Deletions[C_S].NT_size << " \"" << Deletions[C_S].NT_str << "\""
                << "\tChrID " << Deletions[C_S].FragName
                << "\tBP " << Deletions[C_S].BPLeft + 1 << "\t" << Deletions[C_S].BPRight + 1
                << "\tBP_range " << RealStart + 1 << "\t" << RealEnd + 1
                << "\tSupports " << NumberOfReads
                << "\t+ " << LeftS - 1
                << "\t- " << RightS - 1
                << "\tS1 " << EasyScore << "\tS2 " << PreciseScore;

   int SUM_MS = 0;
   for (unsigned int i = C_S; i <= C_E; i++) SUM_MS += Deletions[i].MS;
   DeletionOutf << "\tSUM_MS " << SUM_MS;

   short NumberSupSamples = 0;
   for (short i = 0; i < VectorTag.size(); ++i)
     if (NumSupportPerTag[i]) ++NumberSupSamples;
   DeletionOutf << "\tNumSupSamples " << NumberSupSamples;
   for (short i = 0; i < VectorTag.size(); i++)
     if (NumSupportPerTag[i])
        DeletionOutf << "\t" << VectorTag[i] << " " << NumSupportPerTag[i];
   DeletionOutf << endl;

   //DeletionOutf << TheInput.substr(Deletions[C_S].Left - ReportLength + Deletions[C_S].BP + 1, ReportLength * 2) << endl;// << endl;// ReportLength
   DeletionOutf << TheInput.substr(Deletions[C_S].Left - ReportLength + Deletions[C_S].BP + 1, ReportLength);// << endl;// ReportLength
   if (Deletions[C_S].IndelSize >= 14) {
      DeletionOutf << Cap2Low(TheInput.substr(Deletions[C_S].Left + Deletions[C_S].BP + 1, 5))
                   << "<" << Deletions[C_S].IndelSize - 10 << ">"
                   << Cap2Low(TheInput.substr(Deletions[C_S].Right - Deletions[C_S].ReadLength + Deletions[C_S].BP - 3, 5));
   }
   else DeletionOutf << Cap2Low(TheInput.substr(Deletions[C_S].Left + Deletions[C_S].BP + 1, GapSize));
   DeletionOutf << TheInput.substr(Deletions[C_S].Left + Deletions[C_S].BP + 1 + Deletions[C_S].IndelSize, ReportLength - GapSize) << endl;// ReportLength
   short SpaceBeforeReadSeq;
   for (unsigned int GoodIndex = C_S; GoodIndex <= C_E; GoodIndex++) {
      SpaceBeforeReadSeq = ReportLength - Deletions[GoodIndex].BP - 1;

      for (int i = 0; i < SpaceBeforeReadSeq; i++) DeletionOutf << " ";
      short SpaceBeforeD = ReportLength + ReportLength - SpaceBeforeReadSeq - Deletions[GoodIndex].ReadLength;
      if (Deletions[GoodIndex].MatchedD == Minus) {
         DeletionOutf << Deletions[GoodIndex].UnmatchedSeq.substr(0, Deletions[GoodIndex].BP + 1);// << endl;
         for (int i = 0; i < GapSize; i++) DeletionOutf << " ";
         DeletionOutf << Deletions[GoodIndex].UnmatchedSeq.substr(Deletions[GoodIndex].BP + 1, Deletions[GoodIndex].ReadLength - Deletions[GoodIndex].BP);// << endl;
      }
      else {
         DeletionOutf << ReverseComplement(Deletions[GoodIndex].UnmatchedSeq).substr(0, Deletions[GoodIndex].BP + 1);// << endl;
         for (int i = 0; i < GapSize; i++) DeletionOutf << " ";
         DeletionOutf << ReverseComplement(Deletions[GoodIndex].UnmatchedSeq).substr(Deletions[GoodIndex].BP + 1, Deletions[GoodIndex].ReadLength - Deletions[GoodIndex].BP);// << endl;
      }
      for (int i = 0; i < SpaceBeforeD; i++) DeletionOutf << " ";
      DeletionOutf << "\t" << Deletions[GoodIndex].MatchedD << "\t"
                   << Deletions[GoodIndex].MatchedRelPos
		   << "\t" << Deletions[GoodIndex].MS
                   << "\t" << Deletions[GoodIndex].Tag
		   << "\t" <<  Deletions[GoodIndex].Name << endl;
                   //<< "\t" << Deletions[C_S].BPLeft
                   //<< "\t" << Deletions[C_S].BPRight << endl;
   }
}

void OutputInversions(const vector <SPLIT_READ> & Inv,
							 const string & TheInput,
							 const unsigned int & C_S,
							 const unsigned int & C_E,
							 const unsigned int & RealStart,
							 const unsigned int & RealEnd,
							 ofstream & InvOutf) {
	int LeftNT_index = -1;
	int RightNT_index = -1;
	for (unsigned Index = C_S; Index <= C_E; Index++) {
		if (Inv[Index].MatchedD == Plus) {
			LeftNT_index = Index;
			break;
		}
	}
	for (unsigned Index = C_S; Index <= C_E; Index++) {
		if (Inv[Index].MatchedD == Minus) {
			RightNT_index = Index;
			break;
		}
	}
	short LeftNT_size = 0;
	short RightNT_size = 0;
	string LeftNT_str = "";
	string RightNT_str = "";
	if (LeftNT_index != -1) {
		LeftNT_size = Inv[LeftNT_index].NT_size;
		LeftNT_str = Inv[LeftNT_index].NT_str;
	}
	if (RightNT_index != -1) {
		RightNT_size = Inv[RightNT_index].NT_size;
		RightNT_str = Inv[RightNT_index].NT_str;
	}
   unsigned int NumberOfReads = C_E - C_S + 1;
   float LeftScore = 0;
   float RightScore = 0;
   unsigned int LeftS = 1;
   unsigned int RightS = 1;
   unsigned int NumSupportPerTag[VectorTag.size()];
   for (short i = 0; i < VectorTag.size(); i++) NumSupportPerTag[i] = 0;
   for (unsigned int i = C_S; i <= C_E; i++) {
      for (short j = 0; j < VectorTag.size(); j++) {
			if (Inv[i].Tag == VectorTag[j]) NumSupportPerTag[j]++;
      }
      if (Inv[i].MatchedD == Plus) {
         LeftScore += Inv[i].score;
         LeftS++;
      }
      if (Inv[i].MatchedD == Minus) {
         RightScore += Inv[i].score;
         RightS++;
      }
   }
   unsigned int EasyScore = LeftS * RightS;
   double PreciseScore = (LeftScore + RightScore) * (-1);
   short GapSize =0;
   if (Inv[C_S].IndelSize < 14) GapSize = Inv[C_S].IndelSize;
   else GapSize = 13 + (int)log10(Inv[C_S].IndelSize - 10);

   InvOutf << "####################################################################################################" << endl;
   InvOutf << NumberOfInvInstances << "\tINV " << Inv[C_S].IndelSize// << " bases "
	<< "\tNT " <<  LeftNT_size << ":" << RightNT_size << " \"" << LeftNT_str << "\":\"" << RightNT_str << "\""
	<< "\tChrID " << Inv[C_S].FragName
	<< "\tBP " << Inv[C_S].BPLeft + 1 << "\t" << Inv[C_S].BPRight + 1
	<< "\tBP_range " << Inv[C_S].BPLeft + 1 << "\t" << Inv[C_S].BPRight + 1
	<< "\tSupports " << NumberOfReads
	<< "\t+ " << LeftS - 1
	<< "\t- " << RightS - 1
	<< "\tS1 " << EasyScore << "\tS2 " << PreciseScore;

   int SUM_MS = 0;
   for (unsigned int i = C_S; i <= C_E; i++) SUM_MS += Inv[i].MS;
   InvOutf << "\tSUM_MS " << SUM_MS;

   short NumberSupSamples = 0;
   for (short i = 0; i < VectorTag.size(); ++i)
		if (NumSupportPerTag[i]) ++NumberSupSamples;
   InvOutf << "\tNumSupSamples " << NumberSupSamples;
   for (short i = 0; i < VectorTag.size(); i++)
		if (NumSupportPerTag[i])
			InvOutf << "\t" << VectorTag[i] << " " << NumSupportPerTag[i];
   InvOutf << endl;

	short SpaceBeforeReadSeq;
   //DeletionOutf << TheInput.substr(Deletions[C_S].Left - ReportLength + Deletions[C_S].BP + 1, ReportLength * 2) << endl;// << endl;// ReportLength
   InvOutf << TheInput.substr(Inv[C_S].BPLeft + SpacerBeforeAfter - ReportLength, ReportLength);//;// << endl;// ReportLength
	//cout << Inv[C_S].NT_size << "\t" << Inv[C_S].NT_2size << endl;
	if (LeftNT_size) {
		for (int i = 0; i < LeftNT_size; i++) {
			InvOutf << " ";
		}
	}
	InvOutf << Cap2Low( ReverseComplement(TheInput.substr(Inv[C_S].BPRight + 1 + SpacerBeforeAfter - ReportLength, ReportLength))) << endl;
	for (unsigned int GoodIndex = C_S; GoodIndex <= C_E; GoodIndex++) {
      if (Inv[GoodIndex].MatchedD == Plus && Inv[GoodIndex].MatchedRelPos < Inv[GoodIndex].BPLeft) {
			SpaceBeforeReadSeq = ReportLength - Inv[GoodIndex].BP - 1;
			for (int i = 0; i < SpaceBeforeReadSeq; i++) InvOutf << " ";
         InvOutf << ReverseComplement(Inv[GoodIndex].UnmatchedSeq);
			for (int i = 0; i < Inv[GoodIndex].BP; i++) InvOutf << " ";
			InvOutf //<< "\t" << Inv[GoodIndex].NT_size << "\t\"" << Inv[GoodIndex].NT_str
			        << "\t" << Inv[GoodIndex].MatchedD << "\t"
			        << Inv[GoodIndex].MatchedRelPos
			        << "\t" << Inv[GoodIndex].MS
			        << "\t" << Inv[GoodIndex].Tag
			        << "\t" <<  Inv[GoodIndex].Name << endl;
			//<< "\t" << Deletions[C_S].BPLeft
			//<< "\t" << Deletions[C_S].BPRight << endl;
      }
      else if (Inv[GoodIndex].MatchedD == Plus && Inv[GoodIndex].MatchedRelPos > Inv[GoodIndex].BPLeft) {
			SpaceBeforeReadSeq = ReportLength - Inv[GoodIndex].BP - 1;
			for (int i = 0; i < SpaceBeforeReadSeq; i++) InvOutf << " ";
         InvOutf << Inv[GoodIndex].UnmatchedSeq;
			InvOutf //<< "\t" << Inv[GoodIndex].NT_size << "\t\"" << Inv[GoodIndex].NT_str
			        << "\t" << Inv[GoodIndex].MatchedD << "\t"
			        << Inv[GoodIndex].MatchedRelPos
			        << "\t" << Inv[GoodIndex].MS
			        << "\t" << Inv[GoodIndex].Tag
			        << "\t" <<  Inv[GoodIndex].Name << endl;
      }
   }

	InvOutf << "----------------------------------------------------------------------------------------------------" << endl;


	InvOutf << Cap2Low( ReverseComplement(TheInput.substr(Inv[C_S].BPLeft + SpacerBeforeAfter, ReportLength)));
	if (RightNT_size) {
		for (int i = 0; i < RightNT_size; i++) {
			InvOutf << " ";
		}
	}
	InvOutf << TheInput.substr(Inv[C_S].BPRight + 1 + SpacerBeforeAfter, ReportLength) << endl;
	for (unsigned int GoodIndex = C_S; GoodIndex <= C_E; GoodIndex++) {
      if (Inv[GoodIndex].MatchedD == Minus && Inv[GoodIndex].MatchedRelPos < Inv[GoodIndex].BPRight) {
			SpaceBeforeReadSeq = ReportLength - Inv[GoodIndex].BP - 1;
			for (int i = 0; i < SpaceBeforeReadSeq; i++) InvOutf << " ";
         InvOutf << ReverseComplement(Inv[GoodIndex].UnmatchedSeq);
			for (int i = 0; i < Inv[GoodIndex].BP; i++) InvOutf << " ";
			InvOutf //<< "\t" << Inv[GoodIndex].NT_size << "\t\"" << Inv[GoodIndex].NT_str
			        << "\t" << Inv[GoodIndex].MatchedD << "\t"
			<< Inv[GoodIndex].MatchedRelPos
			<< "\t" << Inv[GoodIndex].MS
			<< "\t" << Inv[GoodIndex].Tag
			<< "\t" <<  Inv[GoodIndex].Name << endl;
			//<< "\t" << Deletions[C_S].BPLeft
			//<< "\t" << Deletions[C_S].BPRight << endl;
      }
      else if (Inv[GoodIndex].MatchedD == Minus && Inv[GoodIndex].MatchedRelPos > Inv[GoodIndex].BPRight) {
			SpaceBeforeReadSeq = ReportLength - Inv[GoodIndex].BP - 1;
			for (int i = 0; i < SpaceBeforeReadSeq; i++) InvOutf << " ";
         InvOutf << Inv[GoodIndex].UnmatchedSeq;
			InvOutf //<< "\t" << Inv[GoodIndex].NT_size << "\t\"" << Inv[GoodIndex].NT_str
			        << "\t" << Inv[GoodIndex].MatchedD << "\t"
			<< Inv[GoodIndex].MatchedRelPos
			<< "\t" << Inv[GoodIndex].MS
			<< "\t" << Inv[GoodIndex].Tag
			<< "\t" <<  Inv[GoodIndex].Name << endl;
      }
   }


}

void OutputSIs(const vector <SPLIT_READ> & SIs,
               const string & TheInput,
               const unsigned int & C_S,
               const unsigned int & C_E,
	       const unsigned int & RealStart,
	       const unsigned int & RealEnd,
               ofstream & SIsOutf) {
   //short ReadLength = SIs[C_S].ReadLength;
   //short ReadLengthMinus = ReadLength - 1;
   unsigned int NumberOfReads = C_E - C_S + 1;
   float LeftScore = 0;
   float RightScore = 0;
   unsigned int LeftS = 1;
   unsigned int RightS = 1;
   unsigned int NumSupportPerTag[VectorTag.size()];
   for (short i = 0; i < VectorTag.size(); i++) NumSupportPerTag[i] = 0;
   for (unsigned int i = C_S; i <= C_E; i++) {
      for (short j = 0; j < VectorTag.size(); j++) {
	if (SIs[i].Tag == VectorTag[j]) NumSupportPerTag[j]++;
      }
      if (SIs[i].MatchedD == Plus) {
         LeftScore += SIs[i].score;
         LeftS++;
      }
      if (SIs[i].MatchedD == Minus) {
         RightScore += SIs[i].score;
         RightS++;
      }
   }
   unsigned int EasyScore = LeftS * RightS;
   double PreciseScore = (LeftScore + RightScore) * (-1);
   string CurrentReadSeq;

   SIsOutf << "####################################################################################################" << endl;
   SIsOutf << NumberOfSIsInstances << "\tI " << SIs[C_S].IndelSize << "\tNT " << SIs[C_S].IndelSize << " \""
           //<< SIs[C_S].InsertedStr
           << GetConsensusInsertedStr(SIs, C_S, C_E)
           << "\""
           << "\tChrID " << SIs[C_S].FragName
           << "\tBP " << SIs[C_S].BPLeft + 1 << "\t" << SIs[C_S].BPRight + 1
           << "\tBP_range " << RealStart + 1 << "\t" << RealEnd + 1
           << "\tSupports " << NumberOfReads
           << "\t+ " << LeftS - 1
           << "\t- " << RightS - 1
           << "\tS1 " << EasyScore << "\tS2 " << PreciseScore;

   int SUM_MS = 0;
   for (unsigned int i = C_S; i <= C_E; i++) SUM_MS += SIs[i].MS;
   SIsOutf << "\tSUM_MS " << SUM_MS;

   short NumberSupSamples = 0;
   for (short i = 0; i < VectorTag.size(); ++i)
     if (NumSupportPerTag[i]) ++NumberSupSamples;
   SIsOutf << "\tNumSupSamples " << NumberSupSamples;
   for (short i = 0; i < VectorTag.size(); i++)
     if (NumSupportPerTag[i])
        SIsOutf << "\t" << VectorTag[i] << " " << NumSupportPerTag[i];
   SIsOutf << endl;

   SIsOutf << TheInput.substr(SIs[C_S].Left - ReportLength + SIs[C_S].BP + 1, ReportLength);// ReportLength
   for (int i = 0; i < SIs[C_S].IndelSize; i++) SIsOutf << " ";
   SIsOutf << TheInput.substr(SIs[C_S].Left + SIs[C_S].BP + 1, ReportLength) << endl;// ReportLength

   for (unsigned int GoodIndex = C_S; GoodIndex <= C_E; GoodIndex++) {
      short SpaceBeforeReadSeq = ReportLength - SIs[GoodIndex].BP - 1;
      for (short i = 0; i < SpaceBeforeReadSeq; i++) SIsOutf << " ";
      if (SIs[GoodIndex].MatchedD == Minus)
         SIsOutf << SIs[GoodIndex].UnmatchedSeq;
      else SIsOutf << ReverseComplement(SIs[GoodIndex].UnmatchedSeq);
      short SpaceBeforeD = ReportLength + ReportLength - SpaceBeforeReadSeq - SIs[GoodIndex].ReadLength;
      for (short i = 0; i < SpaceBeforeD; i++) SIsOutf << " ";
      SIsOutf << "\t" << SIs[GoodIndex].MatchedD
              << "\t" << SIs[GoodIndex].MatchedRelPos
              << "\t" << SIs[GoodIndex].MS
              << "\t" << SIs[GoodIndex].Tag << "\t" << SIs[GoodIndex].Name << endl;
   }
}

void OutputDI(const vector <SPLIT_READ> & DI,
                     const string & TheInput,
                     const unsigned int & C_S,
                     const unsigned int & C_E,
		     const unsigned int & RealStart,
		     const unsigned int & RealEnd,
                     ofstream & DeletionOutf) {
   //short ReadLength = DI[C_S].ReadLength;
   //short ReadLengthMinus = ReadLength - 1;
   unsigned int NumberOfReads = C_E - C_S + 1;
   //float LeftScore = 0;
   //float RightScore = 0;
   unsigned int LeftS = 1;
   unsigned int RightS = 1;
   unsigned int NumSupportPerTag[VectorTag.size()];
   for (short i = 0; i < VectorTag.size(); i++) NumSupportPerTag[i] = 0;
   for (unsigned int i = C_S; i <= C_E; i++) {
      for (short j = 0; j < VectorTag.size(); j++) {
	if (DI[i].Tag == VectorTag[j]) NumSupportPerTag[j]++;
      }
      if (DI[i].MatchedD == Plus) {
         //LeftScore += DI[i].score;
         LeftS++;
      }
      if (DI[i].MatchedD == Minus) {
         //RightScore += DI[i].score;
         RightS++;
      }
   }
   unsigned int EasyScore = LeftS * RightS;
   //double PreciseScore = (LeftScore + RightScore) * (-1);
   //short GapSize =0;
   //if (DI[C_S].IndelSize < 14) GapSize = DI[C_S].IndelSize;
   //else GapSize = 13 + (int)log10(Deletions[C_S].IndelSize - 10);

   DeletionOutf << "####################################################################################################" << endl;
   DeletionOutf << NumberOfDIInstances + NumberOfDeletionsInstances + 1 << "\tD " << DI[C_S].IndelSize
	            << "\tNT " << DI[C_S].NT_size << " \""
                << DI[C_S].NT_str
                //<< GetConsensusInsertedStr(DI, C_S, C_E)
                << "\""//
                << "\tChrID " << DI[C_S].FragName
                << "\tBP " << DI[C_S].BPLeft + 1 << "\t" << DI[C_S].BPRight + 1
                << "\tSupports " << NumberOfReads
                << "\t+ " << LeftS - 1
                << "\t- " << RightS - 1
                << "\tS1 " << EasyScore;// << "\tS2 0.0";// << PreciseScore << "\t";

   int SUM_MS = 0;
   for (unsigned int i = C_S; i <= C_E; i++) SUM_MS += DI[i].MS;
   DeletionOutf << "\tSUM_MS " << SUM_MS;

   short NumberSupSamples = 0;
   for (short i = 0; i < VectorTag.size(); ++i)
     if (NumSupportPerTag[i]) ++NumberSupSamples;
   DeletionOutf << "\tNumSupSamples " << NumberSupSamples;
   for (short i = 0; i < VectorTag.size(); i++)
      if (NumSupportPerTag[i])
         DeletionOutf << "\t" << VectorTag[i] << " " << NumSupportPerTag[i];
   DeletionOutf << endl;

   //DeletionOutf << TheInput.substr(DI[C_S].Left - ReportLength + DI[C_S].BP + 1, 2 * ReportLength) << endl;
   DeletionOutf << TheInput.substr(DI[C_S].Left - ReportLength + DI[C_S].BP + 1, ReportLength);// << endl;// ReportLength

   for (short i = 0; i < DI[C_S].NT_size; i++) DeletionOutf << " ";
   DeletionOutf << TheInput.substr(DI[C_S].Left + DI[C_S].BP + 1 + DI[C_S].IndelSize, ReportLength) << endl;// ReportLength
   short SpaceBeforeReadSeq;
   for (unsigned int GoodIndex = C_S; GoodIndex <= C_E; GoodIndex++) {
      SpaceBeforeReadSeq = ReportLength - DI[GoodIndex].BP - 1;

      for (int i = 0; i < SpaceBeforeReadSeq; i++) DeletionOutf << " ";
      //short SpaceBeforeD = ReportLength + ReportLength - SpaceBeforeReadSeq - Deletions[GoodIndex].ReadLength;
      if (DI[GoodIndex].MatchedD == Minus) {
         DeletionOutf << DI[GoodIndex].UnmatchedSeq << "\t";
         //for (int i = 0; i < GapSize; i++) DeletionOutf << " ";
         //DeletionOutf << Deletions[GoodIndex].UnmatchedSeq.substr(Deletions[GoodIndex].BP + 1, Deletions[GoodIndex].ReadLength - Deletions[GoodIndex].BP);// << endl;
      }
      else {
         DeletionOutf << ReverseComplement(DI[GoodIndex].UnmatchedSeq) << "\t";
         //for (int i = 0; i < GapSize; i++) DeletionOutf << " ";
         //DeletionOutf << ReverseComplement(Deletions[GoodIndex].UnmatchedSeq).substr(Deletions[GoodIndex].BP + 1, Deletions[GoodIndex].ReadLength - Deletions[GoodIndex].BP);// << endl;
      }
      //for (int i = 0; i < SpaceBeforeD; i++) DeletionOutf << " ";
      DeletionOutf << "\t" << DI[GoodIndex].MatchedD << "\t" << DI[GoodIndex].MatchedRelPos
                   << "\t" << DI[GoodIndex].MS
                   << "\t" << DI[GoodIndex].Tag << "\t" <<  DI[GoodIndex].Name << endl;
   }
}


void SortOutputSI(const unsigned & NumBoxes, const string & CurrentChr, vector <SPLIT_READ> & Reads, vector <unsigned> SIs[], ofstream & SIsOutf) {
   cout << "Sorting and outputing short insertions ..." << endl;
   unsigned int SIsNum;
   short CompareResult;
   SPLIT_READ Temp4Exchange;

	vector <SPLIT_READ> InputIndels;
   vector <SPLIT_READ> GoodIndels;
   unsigned int GoodNum;
   vector <Indel4output> IndelEvents;

   for (unsigned Box_index = 0; Box_index < NumBoxes; Box_index++) {
      if (!SIs[Box_index].empty()) {
			InputIndels.clear();
			SIsNum = SIs[Box_index].size();
			//cout << "SIsNum " << SIsNum << endl;
			for (int i = 0; i < SIsNum; i++) {
				InputIndels.push_back(Reads[SIs[Box_index][i]]);
			}
			for (unsigned int First = 0; First < SIsNum - 1; First++) {
            //if (InputIndels[First].OK)
				{
               for (unsigned int Second = First + 1; Second < SIsNum; Second++) {
                  //if (InputIndels[Second].OK)
						{
							if (InputIndels[First].BPLeft < InputIndels[Second].BPLeft) continue;
							else if (InputIndels[First].BPLeft > InputIndels[Second].BPLeft) {
								CompareResult = 1;
							}
							else {
								if (InputIndels[First].IndelSize < InputIndels[Second].IndelSize) continue;
								else if (InputIndels[First].IndelSize > InputIndels[Second].IndelSize) {
									CompareResult = 1;
								}
								else { // InputIndels[First].BPRight == InputIndels[Second].BPRight
										short Compare2Str = CompareTwoString(InputIndels[First].InsertedStr, InputIndels[Second].InsertedStr );
                                    //cout << Compare2Str << " "
                                    //     << First << " " << InputIndels[First].InsertedStr << " "
                                    //     << Second << " " << InputIndels[Second].InsertedStr << endl;
										if (Compare2Str > 0) CompareResult = 1;
										else if (Compare2Str == 0) CompareResult = 2;
										else continue;

								}
								//else if (InputIndels[First].MatchedRelPos == InputIndels[Second].MatchedRelPos) {
								//	if (InputIndels[First].UnmatchedSeq == InputIndels[Second].UnmatchedSeq) {
								//		InputIndels[Second].OK = false;
								//	}
								//}
							}
                     //CompareResult = CompareTwoReads(InputIndels[First], InputIndels[Second]);
                     if (CompareResult == 1) {
  	                     Temp4Exchange = InputIndels[First];
  	                     InputIndels[First] = InputIndels[Second];
  	                     InputIndels[Second] = Temp4Exchange;
                     }
                     //else if (CompareResult == 2) InputIndels[Second].OK = false;
                  }
               }
            }
         }
          //for (unsigned int First = 0; First < SIsNum - 1; First++) {
          //    cout << InputIndels[First].BPLeft << " " << InputIndels[First].IndelSize << " " << InputIndels[First].InsertedStr << endl;
          //}
         GoodIndels.clear();
	      IndelEvents.clear();
			//cout << GoodIndels.size() << endl;
         for (unsigned int First = 0; First < SIsNum; First++) {
            //if (InputIndels[First].OK)
					GoodIndels.push_back(InputIndels[First]);
         }
         GoodNum = GoodIndels.size();
			//cout << Box_index << " " << GoodNum << endl;
			if (GoodNum == 0) continue;
         Indel4output OneIndelEvent;
         OneIndelEvent.Start = 0;
	 OneIndelEvent.End = 0;
         OneIndelEvent.Support = OneIndelEvent.End - OneIndelEvent.Start + 1;
	 OneIndelEvent.IndelSize = GoodIndels[0].IndelSize;
         OneIndelEvent.IndelStr = GoodIndels[0].InsertedStr;
	 OneIndelEvent.BPLeft =  GoodIndels[0].BPLeft;
	 OneIndelEvent.BPRight =  GoodIndels[0].BPRight;
         //OneIndelEvent.IndelStr = GoodIndels[0].InsertedStr;
	     OneIndelEvent.WhetherReport = true;
         for (unsigned int GoodIndex = 1; GoodIndex < GoodNum; GoodIndex++) {
            if (GoodIndels[GoodIndex].BPLeft == OneIndelEvent.BPLeft
                && GoodIndels[GoodIndex].IndelSize == OneIndelEvent.IndelSize
                && OneIndelEvent.IndelStr == GoodIndels[GoodIndex].InsertedStr )

               OneIndelEvent.End = GoodIndex;
	        else  {
                   OneIndelEvent.RealStart = OneIndelEvent.BPLeft;
                   OneIndelEvent.RealEnd = OneIndelEvent.BPRight;

                   OneIndelEvent.Support = OneIndelEvent.End - OneIndelEvent.Start + 1;
                   GetRealStart4Insertion(CurrentChr, OneIndelEvent.IndelStr, OneIndelEvent.RealStart, OneIndelEvent.RealEnd);
	           IndelEvents.push_back(OneIndelEvent);
                   OneIndelEvent.Start = GoodIndex;
		   OneIndelEvent.End = GoodIndex;
		   OneIndelEvent.BPLeft =  GoodIndels[GoodIndex].BPLeft;
                   OneIndelEvent.BPRight =  GoodIndels[GoodIndex].BPRight;
		   OneIndelEvent.IndelSize =  GoodIndels[GoodIndex].IndelSize;
                   OneIndelEvent.IndelStr =  GoodIndels[GoodIndex].InsertedStr;
                   //OneIndelEvent.IndelStr = GoodIndels[GoodIndex].InsertedStr;
            }
         }

	     //if (OneIndelEvent.End - OneIndelEvent.Start + 1 >= NumRead2ReportCutOff)
             OneIndelEvent.RealStart = OneIndelEvent.BPLeft;
             OneIndelEvent.RealEnd = OneIndelEvent.BPRight;

             OneIndelEvent.Support = OneIndelEvent.End - OneIndelEvent.Start + 1;
             GetRealStart4Insertion(CurrentChr, OneIndelEvent.IndelStr, OneIndelEvent.RealStart, OneIndelEvent.RealEnd);
	     IndelEvents.push_back(OneIndelEvent);

	     //string IndelType;
	     unsigned int RealStart;
	     unsigned int RealEnd;
	     //bool WhetherDeletion = true;
	     string IndelStr;
	     unsigned int Max_Support;
	     unsigned int Max_Support_Index;
             unsigned int IndelSize;
	    if (IndelEvents.size()) {
            for (unsigned EventIndex = 0; EventIndex < IndelEvents.size(); EventIndex++) {
               if (IndelEvents[EventIndex].WhetherReport) {
                  RealStart = IndelEvents[EventIndex].RealStart;
                  RealEnd = IndelEvents[EventIndex].RealEnd;
                  IndelSize = IndelEvents[EventIndex].IndelSize;
                  Max_Support = IndelEvents[EventIndex].Support;
                  Max_Support_Index = EventIndex;

                  for (unsigned EventIndex_left = 0; EventIndex_left < IndelEvents.size(); EventIndex_left++) {
                     if (IndelEvents[EventIndex_left].WhetherReport == false) continue;
                     else if (IndelEvents[EventIndex_left].RealStart != RealStart) continue;
                     else if (IndelEvents[EventIndex_left].RealEnd != RealEnd) continue;
                     else if (IndelEvents[EventIndex_left].IndelSize != IndelSize) continue;
                     else {
                        IndelEvents[EventIndex_left].WhetherReport = false;
                        if (IndelEvents[EventIndex_left].Support > Max_Support) {
                           Max_Support = IndelEvents[EventIndex_left].Support;
                           Max_Support_Index = EventIndex_left;
                        }
                     }
                  }
                  // report max one
                  //cout << Max_Support << endl;
                  if (Max_Support >= NumRead2ReportCutOff) {
		     OutputSIs(GoodIndels, CurrentChr,
		               IndelEvents[Max_Support_Index].Start,
			       IndelEvents[Max_Support_Index].End,
			       RealStart, RealEnd, SIsOutf);
		     NumberOfSIsInstances++;
		  }
               }
            }
         }
      }   // if (!insertion[Box_index].empty())
   }
   cout << "Short insertions: " << NumberOfSIsInstances << endl << endl;
}

void SortOutputTD(const unsigned & NumBoxes, const string & CurrentChr, vector <SPLIT_READ> & AllReads, vector <unsigned> TDs[], ofstream & TDOutf) {
   cout << "Sorting and outputing tandem duplications ..." << endl;
   unsigned int TDNum;
   short CompareResult;
   SPLIT_READ Temp4Exchange;

   unsigned int GoodNum;
	vector <SPLIT_READ> InputIndels;
   vector <SPLIT_READ> GoodIndels;
   vector <Indel4output> IndelEvents;

   for (unsigned Box_index = 0; Box_index < NumBoxes; Box_index++) {
      if (!TDs[Box_index].empty()) {
			InputIndels.clear();
			TDNum = TDs[Box_index].size();
			for (int i = 0 ; i < TDNum; i ++) {
				InputIndels.push_back(AllReads[TDs[Box_index][i]]);
			}

         for (unsigned int First = 0; First < TDNum - 1; First++) {
            //if (InputIndels[First].OK)
				{
               for (unsigned int Second = First + 1; Second < TDNum; Second++) {
                  //if (InputIndels[Second].OK)
						{
							if (InputIndels[First].BPLeft < InputIndels[Second].BPLeft) continue;
                     else if (InputIndels[First].BPLeft > InputIndels[Second].BPLeft) {
								CompareResult = 1;
							}
							else if (InputIndels[First].BPLeft == InputIndels[Second].BPLeft) {
								if (InputIndels[First].BPRight < InputIndels[Second].BPRight) continue;
								else if (InputIndels[First].BPRight > InputIndels[Second].BPRight) {
									CompareResult = 1;
								}
								//else {
								//	if (InputIndels[First].MatchedRelPos == InputIndels[Second].MatchedRelPos) {
								//		if (InputIndels[First].UnmatchedSeq == InputIndels[Second].UnmatchedSeq) {
								//			InputIndels[Second].OK = false;
								//		}
								//
								//	}
								//}
							}
							if (CompareResult == 1) {
								Temp4Exchange = InputIndels[First];
  	                     InputIndels[First] = InputIndels[Second];
  	                     InputIndels[Second] = Temp4Exchange;
							}
                  }
               }
            }
         }
         GoodIndels.clear();
			IndelEvents.clear();

         for (unsigned int First = 0; First < TDNum; First++) {
            //if (InputIndels[First].OK)
					GoodIndels.push_back(InputIndels[First]);
         }
         GoodNum = GoodIndels.size();
			//cout << Box_index << " " << GoodNum << endl;
			if (GoodNum == 0) continue;
			//    cout << GoodNum << endl;
         Indel4output OneIndelEvent;
         OneIndelEvent.Start = 0;
			OneIndelEvent.End = 0;
         OneIndelEvent.Support = OneIndelEvent.End - OneIndelEvent.Start + 1;
			OneIndelEvent.BPLeft =  GoodIndels[0].BPLeft;
         OneIndelEvent.BPRight =  GoodIndels[0].BPRight;
			OneIndelEvent.WhetherReport = true;
         for (unsigned int GoodIndex = 1; GoodIndex < GoodNum; GoodIndex++) {
            if (GoodIndels[GoodIndex].BPLeft == OneIndelEvent.BPLeft
                && GoodIndels[GoodIndex].BPRight == OneIndelEvent.BPRight)
               OneIndelEvent.End = GoodIndex;
				else  {
					OneIndelEvent.RealStart = OneIndelEvent.BPLeft;
					OneIndelEvent.RealEnd = OneIndelEvent.BPRight;
					OneIndelEvent.Support = OneIndelEvent.End - OneIndelEvent.Start + 1;
					GetRealStart4Deletion(CurrentChr, OneIndelEvent.RealStart, OneIndelEvent.RealEnd);
					IndelEvents.push_back(OneIndelEvent);
					OneIndelEvent.Start = GoodIndex;
					OneIndelEvent.End = GoodIndex;
					OneIndelEvent.BPLeft =  GoodIndels[GoodIndex].BPLeft;
					OneIndelEvent.BPRight =  GoodIndels[GoodIndex].BPRight;
            }
         }

         OneIndelEvent.RealStart = OneIndelEvent.BPLeft;
         OneIndelEvent.RealEnd = OneIndelEvent.BPRight;
         OneIndelEvent.Support = OneIndelEvent.End - OneIndelEvent.Start + 1;
         GetRealStart4Deletion(CurrentChr, OneIndelEvent.RealStart, OneIndelEvent.RealEnd);
         IndelEvents.push_back(OneIndelEvent);
			//	     cout << "IndelEvent: " << IndelEvents.size() << endl;
			//string IndelType;
			unsigned int RealStart;
			unsigned int RealEnd;
			//bool WhetherDeletion = true;
			string IndelStr;
			unsigned int Max_Support;
			unsigned int Max_Support_Index;

			if (IndelEvents.size()) {
            for (unsigned EventIndex = 0; EventIndex < IndelEvents.size(); EventIndex++) {
               if (IndelEvents[EventIndex].WhetherReport) {
                  RealStart = IndelEvents[EventIndex].RealStart;
                  RealEnd = IndelEvents[EventIndex].RealEnd;
                  Max_Support = IndelEvents[EventIndex].Support;
                  Max_Support_Index = EventIndex;
                  for (unsigned EventIndex_left = 0; EventIndex_left < IndelEvents.size(); EventIndex_left++) {
                     if (IndelEvents[EventIndex_left].WhetherReport == false) continue;
                     else if (IndelEvents[EventIndex_left].RealStart != RealStart) continue;
                     else if (IndelEvents[EventIndex_left].RealEnd != RealEnd) continue;
                     else {
                        IndelEvents[EventIndex_left].WhetherReport = false;
                        if (IndelEvents[EventIndex_left].Support > Max_Support) {
                           Max_Support = IndelEvents[EventIndex_left].Support;
                           Max_Support_Index = EventIndex_left;
                        }
                     }
                  }
                  // report max one
                  if (IndelEvents[Max_Support_Index].Support >= NumRead2ReportCutOff) {
							if (GoodIndels[IndelEvents[Max_Support_Index].Start].IndelSize < BalanceCutoff) {
								OutputTDs(GoodIndels, CurrentChr,
													 IndelEvents[Max_Support_Index].Start,
													 IndelEvents[Max_Support_Index].End,
													 RealStart, RealEnd, TDOutf);
								NumberOfTDInstances++;
							}
							else if (ReportEvent(GoodIndels, IndelEvents[Max_Support_Index].Start, IndelEvents[Max_Support_Index].End)) {
								OutputTDs(GoodIndels, CurrentChr,
													 IndelEvents[Max_Support_Index].Start,
													 IndelEvents[Max_Support_Index].End,
													 RealStart, RealEnd, TDOutf);
								NumberOfTDInstances++;
							}
						}
               }
            }
         }
      }   // if (!Deletions[Box_index].empty())
   } // for (unsigned Box_index = 0; Box_index < NumBoxes; Box_index++)
   cout << "Tandem duplications: " << NumberOfTDInstances << endl << endl;
}

void SortOutputTD_NT(const unsigned & NumBoxes, const string & CurrentChr, vector <SPLIT_READ> & AllReads, vector <unsigned> TDs[], ofstream & TDOutf) {
   cout << "Sorting and outputing tandem duplications with non-template sequence ..." << endl;
   unsigned int TDNum;
   short CompareResult;
   SPLIT_READ Temp4Exchange;

	int Count_TD_NT_output = 0;

   unsigned int GoodNum;
	vector <SPLIT_READ> InputIndels;
   vector <SPLIT_READ> GoodIndels;
   vector <Indel4output> IndelEvents;

   for (unsigned Box_index = 0; Box_index < NumBoxes; Box_index++) {
      if (!TDs[Box_index].empty()) {
			InputIndels.clear();
			TDNum = TDs[Box_index].size();
			for (int i = 0 ; i < TDNum; i ++) {
				InputIndels.push_back(AllReads[TDs[Box_index][i]]);
			}

         for (unsigned int First = 0; First < TDNum - 1; First++) {
            //if (InputIndels[First].OK)
				{
               for (unsigned int Second = First + 1; Second < TDNum; Second++) {
                  //if (InputIndels[Second].OK)
						{
							if (InputIndels[First].BPLeft < InputIndels[Second].BPLeft) continue;
                     else if (InputIndels[First].BPLeft > InputIndels[Second].BPLeft) {
								CompareResult = 1;
							}
							else if (InputIndels[First].BPLeft == InputIndels[Second].BPLeft) {
								if (InputIndels[First].BPRight < InputIndels[Second].BPRight) continue;
								else if (InputIndels[First].BPRight > InputIndels[Second].BPRight) {
									CompareResult = 1;
								}
								else { // InputIndels[First].BPRight == InputIndels[Second].BPRight
									if (InputIndels[First].NT_size < InputIndels[Second].NT_size) continue;
									else if (InputIndels[First].NT_size > InputIndels[Second].NT_size) CompareResult = 1;
									else { // InputIndels[First].NT_size == InputIndels[Second].NT_size
										short Compare2Str = CompareTwoString(InputIndels[First].NT_str, InputIndels[Second].NT_str );
										if (Compare2Str > 0) CompareResult = 1;
										else if (Compare2Str == 0) CompareResult = 2;
										else continue;
									}
								}
							}
							if (CompareResult == 1) {
								Temp4Exchange = InputIndels[First];
  	                     InputIndels[First] = InputIndels[Second];
  	                     InputIndels[Second] = Temp4Exchange;
							}
							else if (CompareResult == 2) {
								Temp4Exchange = InputIndels[First + 1];
  	                     InputIndels[First + 1] = InputIndels[Second];
  	                     InputIndels[Second] = Temp4Exchange;
							}
                  }
               }
            }
         }
         GoodIndels.clear();
			IndelEvents.clear();

         for (unsigned int First = 0; First < TDNum; First++) {
            //if (InputIndels[First].OK)
					GoodIndels.push_back(InputIndels[First]);
         }
         GoodNum = GoodIndels.size();
			//cout << Box_index << " " << GoodNum << endl;
			if (GoodNum == 0) continue;
			//    cout << GoodNum << endl;
         Indel4output OneIndelEvent;
         OneIndelEvent.Start = 0;
			OneIndelEvent.End = 0;
         OneIndelEvent.Support = OneIndelEvent.End - OneIndelEvent.Start + 1;
			OneIndelEvent.BPLeft =  GoodIndels[0].BPLeft;
         OneIndelEvent.BPRight =  GoodIndels[0].BPRight;
			OneIndelEvent.WhetherReport = true;
         for (unsigned int GoodIndex = 1; GoodIndex < GoodNum; GoodIndex++) {
            if (GoodIndels[GoodIndex].BPLeft == OneIndelEvent.BPLeft
                && GoodIndels[GoodIndex].BPRight == OneIndelEvent.BPRight)
               OneIndelEvent.End = GoodIndex;
				else  {
					OneIndelEvent.RealStart = OneIndelEvent.BPLeft;
					OneIndelEvent.RealEnd = OneIndelEvent.BPRight;
					OneIndelEvent.Support = OneIndelEvent.End - OneIndelEvent.Start + 1;
					GetRealStart4Deletion(CurrentChr, OneIndelEvent.RealStart, OneIndelEvent.RealEnd);
					IndelEvents.push_back(OneIndelEvent);
					OneIndelEvent.Start = GoodIndex;
					OneIndelEvent.End = GoodIndex;
					OneIndelEvent.BPLeft =  GoodIndels[GoodIndex].BPLeft;
					OneIndelEvent.BPRight =  GoodIndels[GoodIndex].BPRight;
            }
         }

         OneIndelEvent.RealStart = OneIndelEvent.BPLeft;
         OneIndelEvent.RealEnd = OneIndelEvent.BPRight;
         OneIndelEvent.Support = OneIndelEvent.End - OneIndelEvent.Start + 1;
         GetRealStart4Deletion(CurrentChr, OneIndelEvent.RealStart, OneIndelEvent.RealEnd);
         IndelEvents.push_back(OneIndelEvent);
			//	     cout << "IndelEvent: " << IndelEvents.size() << endl;
			//string IndelType;
			unsigned int RealStart;
			unsigned int RealEnd;
			//bool WhetherDeletion = true;
			string IndelStr;
			unsigned int Max_Support;
			unsigned int Max_Support_Index;

			if (IndelEvents.size()) {
            for (unsigned EventIndex = 0; EventIndex < IndelEvents.size(); EventIndex++) {
               if (IndelEvents[EventIndex].WhetherReport) {
                  RealStart = IndelEvents[EventIndex].RealStart;
                  RealEnd = IndelEvents[EventIndex].RealEnd;
                  Max_Support = IndelEvents[EventIndex].Support;
                  Max_Support_Index = EventIndex;
                  for (unsigned EventIndex_left = 0; EventIndex_left < IndelEvents.size(); EventIndex_left++) {
                     if (IndelEvents[EventIndex_left].WhetherReport == false) continue;
                     else if (IndelEvents[EventIndex_left].RealStart != RealStart) continue;
                     else if (IndelEvents[EventIndex_left].RealEnd != RealEnd) continue;
                     else {
                        IndelEvents[EventIndex_left].WhetherReport = false;
                        if (IndelEvents[EventIndex_left].Support > Max_Support) {
                           Max_Support = IndelEvents[EventIndex_left].Support;
                           Max_Support_Index = EventIndex_left;
                        }
                     }
                  }
                  // report max one
                  if (IndelEvents[Max_Support_Index].Support >= NumRead2ReportCutOff) {
							if (GoodIndels[IndelEvents[Max_Support_Index].Start].IndelSize < BalanceCutoff) {
								OutputTDs(GoodIndels, CurrentChr,
											 IndelEvents[Max_Support_Index].Start,
											 IndelEvents[Max_Support_Index].End,
											 RealStart, RealEnd, TDOutf);
								NumberOfTDInstances++;
								Count_TD_NT_output++;
							}
							else if (ReportEvent(GoodIndels, IndelEvents[Max_Support_Index].Start, IndelEvents[Max_Support_Index].End)) {
								OutputTDs(GoodIndels, CurrentChr,
											 IndelEvents[Max_Support_Index].Start,
											 IndelEvents[Max_Support_Index].End,
											 RealStart, RealEnd, TDOutf);
								NumberOfTDInstances++;
								Count_TD_NT_output++;
							}
						}
               }
            }
         }
      }   // if (!Deletions[Box_index].empty())
   } // for (unsigned Box_index = 0; Box_index < NumBoxes; Box_index++)
   cout << "Tandem duplications with non-template sequence (TD_NT): " << Count_TD_NT_output << "\n\n";
}

void SortOutputD(const unsigned & NumBoxes, const string & CurrentChr, vector <SPLIT_READ> & Reads, vector <unsigned> Deletions[], ofstream & DeletionOutf) {
   cout << "Sorting and outputing deletions ..." << endl;
   unsigned int DeletionsNum;
   short CompareResult;
   SPLIT_READ Temp4Exchange;

   unsigned int GoodNum;
	vector <SPLIT_READ> InputIndels;
   vector <SPLIT_READ> GoodIndels;
   vector <Indel4output> IndelEvents;

   for (unsigned Box_index = 0; Box_index < NumBoxes; Box_index++) {
      if (!Deletions[Box_index].empty()) {
			InputIndels.clear();
			DeletionsNum = Deletions[Box_index].size();
			for (int i = 0 ; i < DeletionsNum; i ++) {
				InputIndels.push_back(Reads[Deletions[Box_index][i]]);
			}

         for (unsigned int First = 0; First < DeletionsNum - 1; First++) {
            //if (InputIndels[First].OK)
				{
               for (unsigned int Second = First + 1; Second < DeletionsNum; Second++) {
                  //if (InputIndels[Second].OK)
						{
							if (InputIndels[First].BPLeft < InputIndels[Second].BPLeft) continue;
                     else if (InputIndels[First].BPLeft > InputIndels[Second].BPLeft) {
								CompareResult = 1;
							}
							else if (InputIndels[First].BPLeft == InputIndels[Second].BPLeft) {
								if (InputIndels[First].BPRight < InputIndels[Second].BPRight) continue;
								else if (InputIndels[First].BPRight > InputIndels[Second].BPRight) {
									CompareResult = 1;
								}
								else CompareResult = 2;
								//else {
								//	if (InputIndels[First].MatchedRelPos == InputIndels[Second].MatchedRelPos) {
								//		if (InputIndels[First].UnmatchedSeq == InputIndels[Second].UnmatchedSeq) {
								//			InputIndels[Second].OK = false;
								//		}
								//
								//	}
								//}
							}
							if (CompareResult == 1) {
								Temp4Exchange = InputIndels[First];
  	                     InputIndels[First] = InputIndels[Second];
  	                     InputIndels[Second] = Temp4Exchange;
							}
							else if (CompareResult == 2) {
								Temp4Exchange = InputIndels[First + 1];
  	                     InputIndels[First + 1] = InputIndels[Second];
  	                     InputIndels[Second] = Temp4Exchange;
							}
                  }
               }
            }
         }
         GoodIndels.clear();
	     IndelEvents.clear();

         for (unsigned int First = 0; First < DeletionsNum; First++) {
            //if (InputIndels[First].OK)
					GoodIndels.push_back(InputIndels[First]);
         }
         GoodNum = GoodIndels.size();
			//cout << Box_index << " " << GoodNum << endl;
         if (GoodNum == 0) continue;
   //    cout << GoodNum << endl;
         Indel4output OneIndelEvent;
         OneIndelEvent.Start = 0;
	 OneIndelEvent.End = 0;
         OneIndelEvent.Support = OneIndelEvent.End - OneIndelEvent.Start + 1;
	 OneIndelEvent.BPLeft =  GoodIndels[0].BPLeft;
         OneIndelEvent.BPRight =  GoodIndels[0].BPRight;
	 OneIndelEvent.WhetherReport = true;
         for (unsigned int GoodIndex = 1; GoodIndex < GoodNum; GoodIndex++) {
            if (GoodIndels[GoodIndex].BPLeft == OneIndelEvent.BPLeft
                && GoodIndels[GoodIndex].BPRight == OneIndelEvent.BPRight)
               OneIndelEvent.End = GoodIndex;
	    else  {
                   OneIndelEvent.RealStart = OneIndelEvent.BPLeft;
                   OneIndelEvent.RealEnd = OneIndelEvent.BPRight;
                   OneIndelEvent.Support = OneIndelEvent.End - OneIndelEvent.Start + 1;
                   GetRealStart4Deletion(CurrentChr, OneIndelEvent.RealStart, OneIndelEvent.RealEnd);
	           IndelEvents.push_back(OneIndelEvent);
                   OneIndelEvent.Start = GoodIndex;
		   OneIndelEvent.End = GoodIndex;
		   OneIndelEvent.BPLeft =  GoodIndels[GoodIndex].BPLeft;
		   OneIndelEvent.BPRight =  GoodIndels[GoodIndex].BPRight;
            }
         }

         OneIndelEvent.RealStart = OneIndelEvent.BPLeft;
         OneIndelEvent.RealEnd = OneIndelEvent.BPRight;
         OneIndelEvent.Support = OneIndelEvent.End - OneIndelEvent.Start + 1;
         GetRealStart4Deletion(CurrentChr, OneIndelEvent.RealStart, OneIndelEvent.RealEnd);
         IndelEvents.push_back(OneIndelEvent);
//	     cout << "IndelEvent: " << IndelEvents.size() << endl;
	     //string IndelType;
	     unsigned int RealStart;
	     unsigned int RealEnd;
	     //bool WhetherDeletion = true;
	     string IndelStr;
	     unsigned int Max_Support;
	     unsigned int Max_Support_Index;

	    if (IndelEvents.size()) {
            for (unsigned EventIndex = 0; EventIndex < IndelEvents.size(); EventIndex++) {
               if (IndelEvents[EventIndex].WhetherReport) {
                  RealStart = IndelEvents[EventIndex].RealStart;
                  RealEnd = IndelEvents[EventIndex].RealEnd;
                  Max_Support = IndelEvents[EventIndex].Support;
                  Max_Support_Index = EventIndex;
                  for (unsigned EventIndex_left = 0; EventIndex_left < IndelEvents.size(); EventIndex_left++) {
                     if (IndelEvents[EventIndex_left].WhetherReport == false) continue;
                     else if (IndelEvents[EventIndex_left].RealStart != RealStart) continue;
                     else if (IndelEvents[EventIndex_left].RealEnd != RealEnd) continue;
                     else {
                        IndelEvents[EventIndex_left].WhetherReport = false;
                        if (IndelEvents[EventIndex_left].Support > Max_Support) {
                           Max_Support = IndelEvents[EventIndex_left].Support;
                           Max_Support_Index = EventIndex_left;
                        }
                     }
                  }
                  // report max one
                  if (IndelEvents[Max_Support_Index].Support >= NumRead2ReportCutOff) {
		     if (GoodIndels[IndelEvents[Max_Support_Index].Start].IndelSize < BalanceCutoff) {
		             OutputDeletions(GoodIndels, CurrentChr,
		                             IndelEvents[Max_Support_Index].Start,
				             IndelEvents[Max_Support_Index].End,
				             RealStart, RealEnd, DeletionOutf);
		             NumberOfDeletionsInstances++;
	             }
		     else if (ReportEvent(GoodIndels, IndelEvents[Max_Support_Index].Start, IndelEvents[Max_Support_Index].End)) {
		             OutputDeletions(GoodIndels, CurrentChr,
		                             IndelEvents[Max_Support_Index].Start,
				             IndelEvents[Max_Support_Index].End,
				             RealStart, RealEnd, DeletionOutf);
		             NumberOfDeletionsInstances++;
		     }
		  }
               }
            }
         }
      }   // if (!Deletions[Box_index].empty())
   } // for (unsigned Box_index = 0; Box_index < NumBoxes; Box_index++)
   cout << "Deletions: " << NumberOfDeletionsInstances << endl << endl;
}

void SortOutputInv(const unsigned & NumBoxes, const string & CurrentChr, vector <SPLIT_READ> & Reads, vector <unsigned> Inv[], ofstream & InvOutf) {
   cout << "Sorting and outputing Inversions ..." << endl;
   unsigned int InversionsNum;
   short CompareResult;
   SPLIT_READ Temp4Exchange;
	/*
   unsigned int C_S = 0;
   unsigned int C_E = 0;
   unsigned int C_BP_Left;// = GoodSIs[0].BPLeft;
   unsigned int C_BP_Right;// = GoodSIs[0].BPRight;
   unsigned int C_Indelsize;
	 */
   unsigned int GoodNum;
	vector <SPLIT_READ> InputIndels;
   vector <SPLIT_READ> GoodIndels;
   vector <Indel4output> IndelEvents;

   for (unsigned Box_index = 0; Box_index < NumBoxes; Box_index++) {

      if (!Inv[Box_index].empty()) {
			//cout << Box_index << "\t" << Inv[Box_index].size() << endl;
         InversionsNum = Inv[Box_index].size();
			InputIndels.clear();
			for (int i = 0; i < InversionsNum; i++) InputIndels.push_back(Reads[Inv[Box_index][i]]);
         for (unsigned int First = 0; First < InversionsNum - 1; First++) {
            //if (InputIndels[First].OK)
				{
               for (unsigned int Second = First + 1; Second < InversionsNum; Second++) {
                  //if (InputIndels[Second].OK)
						{
							if (InputIndels[First].BPLeft < InputIndels[Second].BPLeft) continue;
                     else if (InputIndels[First].BPLeft > InputIndels[Second].BPLeft) {
								CompareResult = 1;
							}
							else if (InputIndels[First].BPLeft == InputIndels[Second].BPLeft) {
								if (InputIndels[First].BPRight < InputIndels[Second].BPRight) continue;
								else if (InputIndels[First].BPRight > InputIndels[Second].BPRight) {
									CompareResult = 1;
								}
								else CompareResult = 2;
								//else {
								//	if (InputIndels[First].MatchedRelPos == InputIndels[Second].MatchedRelPos) {
								//		if (InputIndels[First].UnmatchedSeq == InputIndels[Second].UnmatchedSeq) {
								//			InputIndels[Second].OK = false;
								//		}
								//
								//	}
								//}
							}
							if (CompareResult == 1) {
								Temp4Exchange = InputIndels[First];
  	                     InputIndels[First] = InputIndels[Second];
  	                     InputIndels[Second] = Temp4Exchange;
							}
							else if (CompareResult == 2) {
								Temp4Exchange = InputIndels[First + 1];
  	                     InputIndels[First + 1] = InputIndels[Second];
  	                     InputIndels[Second] = Temp4Exchange;
							}
                  }
               }
            }
         }
         GoodIndels.clear();
			IndelEvents.clear();

         for (unsigned int First = 0; First < InversionsNum; First++) {
            //if (InputIndels[First].OK)
					GoodIndels.push_back(InputIndels[First]);
         }
         GoodNum = GoodIndels.size();
			//cout << Box_index << " " << GoodNum << endl;
			if (GoodNum == 0) continue;
			//    cout << GoodNum << endl;
         Indel4output OneIndelEvent;
         OneIndelEvent.Start = 0;
			OneIndelEvent.End = 0;
         OneIndelEvent.Support = OneIndelEvent.End - OneIndelEvent.Start + 1;
			OneIndelEvent.BPLeft =  GoodIndels[0].BPLeft;
         OneIndelEvent.BPRight =  GoodIndels[0].BPRight;
			OneIndelEvent.WhetherReport = true;
         for (unsigned int GoodIndex = 1; GoodIndex < GoodNum; GoodIndex++) {
            if (GoodIndels[GoodIndex].BPLeft == OneIndelEvent.BPLeft
                && GoodIndels[GoodIndex].BPRight == OneIndelEvent.BPRight)
               OneIndelEvent.End = GoodIndex;
				else  {
					OneIndelEvent.RealStart = OneIndelEvent.BPLeft;
					OneIndelEvent.RealEnd = OneIndelEvent.BPRight;
					OneIndelEvent.Support = OneIndelEvent.End - OneIndelEvent.Start + 1;
					//GetRealStart4Deletion(CurrentChr, OneIndelEvent.RealStart, OneIndelEvent.RealEnd);
					IndelEvents.push_back(OneIndelEvent);
					OneIndelEvent.Start = GoodIndex;
					OneIndelEvent.End = GoodIndex;
					OneIndelEvent.BPLeft =  GoodIndels[GoodIndex].BPLeft;
					OneIndelEvent.BPRight =  GoodIndels[GoodIndex].BPRight;
            }
         }

         OneIndelEvent.RealStart = OneIndelEvent.BPLeft;
         OneIndelEvent.RealEnd = OneIndelEvent.BPRight;
         OneIndelEvent.Support = OneIndelEvent.End - OneIndelEvent.Start + 1;
         //GetRealStart4Deletion(CurrentChr, OneIndelEvent.RealStart, OneIndelEvent.RealEnd);
         IndelEvents.push_back(OneIndelEvent);
			//	     cout << "IndelEvent: " << IndelEvents.size() << endl;
			//string IndelType;
			unsigned int RealStart;
			unsigned int RealEnd;
			//bool WhetherDeletion = true;
			string IndelStr;
			unsigned int Max_Support;
			unsigned int Max_Support_Index;

			if (IndelEvents.size()) {
            for (unsigned EventIndex = 0; EventIndex < IndelEvents.size(); EventIndex++) {
					RealStart = IndelEvents[EventIndex].RealStart;
					RealEnd = IndelEvents[EventIndex].RealEnd;
					if (IndelEvents[EventIndex].Support < NumRead2ReportCutOff) continue;
					// report max one
					if (GoodIndels[IndelEvents[EventIndex].Start].IndelSize < BalanceCutoff) {
						OutputInversions(GoodIndels, CurrentChr,
											  IndelEvents[EventIndex].Start,
											  IndelEvents[EventIndex].End,
											  RealStart, RealEnd, InvOutf);
						NumberOfInvInstances++;
					}
					else if (ReportEvent(GoodIndels, IndelEvents[EventIndex].Start, IndelEvents[EventIndex].End)) {
						OutputInversions(GoodIndels, CurrentChr,
											 IndelEvents[EventIndex].Start,
											 IndelEvents[EventIndex].End,
											 RealStart, RealEnd, InvOutf);
						NumberOfInvInstances++;
					}
            }
         }
      }   // if (!Deletions[Box_index].empty())
   } // for (unsigned Box_index = 0; Box_index < NumBoxes; Box_index++)
   cout << "Inversions (INV): " << NumberOfInvInstances << endl << endl;
}

void SortOutputInv_NT(const unsigned & NumBoxes, const string & CurrentChr, vector <SPLIT_READ> & Reads, vector <unsigned> Inv[], ofstream & InvOutf) {
   cout << "Sorting and outputing Inversions with non-template sequence ..." << endl;
   unsigned int InversionsNum;
   short CompareResult;
   SPLIT_READ Temp4Exchange;

	int Count_INV_NT_output = 0;
	/*
	 unsigned int C_S = 0;
	 unsigned int C_E = 0;
	 unsigned int C_BP_Left;// = GoodSIs[0].BPLeft;
	 unsigned int C_BP_Right;// = GoodSIs[0].BPRight;
	 unsigned int C_Indelsize;
	 */
   unsigned int GoodNum;
	vector <SPLIT_READ> InputIndels;
   vector <SPLIT_READ> GoodIndels;
   vector <Indel4output> IndelEvents;

   for (unsigned Box_index = 0; Box_index < NumBoxes; Box_index++) {

      if (!Inv[Box_index].empty()) {
         InversionsNum = Inv[Box_index].size();
			//cout << Box_index << "\t" << Inv[Box_index].size() << endl;
			InputIndels.clear();
			for (int i = 0; i < InversionsNum; i++) InputIndels.push_back(Reads[Inv[Box_index][i]]);
         for (unsigned int First = 0; First < InversionsNum - 1; First++) {

            //if (InputIndels[First].OK)
				{
               for (unsigned int Second = First + 1; Second < InversionsNum; Second++) {
						//cout << InputIndels[First].BPLeft << "\t" << InputIndels[First].BPRight << "\t"
						//     << InputIndels[Second].BPLeft << "\t" << InputIndels[Second].BPRight<< endl;
                  //if (InputIndels[Second].OK)
						{
							if (InputIndels[First].BPLeft < InputIndels[Second].BPLeft) continue;
                     else if (InputIndels[First].BPLeft > InputIndels[Second].BPLeft) {
								CompareResult = 1;
							}
							else if (InputIndels[First].BPLeft == InputIndels[Second].BPLeft) {
								if (InputIndels[First].BPRight < InputIndels[Second].BPRight) continue;
								else if (InputIndels[First].BPRight > InputIndels[Second].BPRight) {
									CompareResult = 1;
								}
								else { // InputIndels[First].BPRight == InputIndels[Second].BPRight
									if (InputIndels[First].NT_size < InputIndels[Second].NT_size) continue;
									else if (InputIndels[First].NT_size > InputIndels[Second].NT_size) CompareResult = 1;
									else { // InputIndels[First].NT_size == InputIndels[Second].NT_size
										short Compare2Str = CompareTwoString(InputIndels[First].NT_str, InputIndels[Second].NT_str );
										if (Compare2Str > 0) CompareResult = 1;
										else if (Compare2Str == 0) CompareResult = 2;
										else continue;
									}
								}
								//else {
								//	if (InputIndels[First].MatchedRelPos == InputIndels[Second].MatchedRelPos) {
								//		if (InputIndels[First].UnmatchedSeq == InputIndels[Second].UnmatchedSeq) {
								//			InputIndels[Second].OK = false;
								//		}
								//
								//	}
								//}
							}
							if (CompareResult == 1) {
								Temp4Exchange = InputIndels[First];
  	                     InputIndels[First] = InputIndels[Second];
  	                     InputIndels[Second] = Temp4Exchange;
							}
							if (CompareResult == 2) {
								Temp4Exchange = InputIndels[First + 1];
  	                     InputIndels[First + 1] = InputIndels[Second];
  	                     InputIndels[Second] = Temp4Exchange;
							}
                  }
               }
            }
         }
         GoodIndels.clear();
			IndelEvents.clear();

         for (unsigned int First = 0; First < InversionsNum; First++) {
            //if (InputIndels[First].OK)
					GoodIndels.push_back(InputIndels[First]);
         }
         GoodNum = GoodIndels.size();
			//cout << "GoodNum " << Box_index << " " << GoodNum << endl;
			if (GoodNum == 0) continue;
			//    cout << GoodNum << endl;
         Indel4output OneIndelEvent;
         OneIndelEvent.Start = 0;
			OneIndelEvent.End = 0;
         OneIndelEvent.Support = OneIndelEvent.End - OneIndelEvent.Start + 1;
			OneIndelEvent.BPLeft =  GoodIndels[0].BPLeft;
         OneIndelEvent.BPRight =  GoodIndels[0].BPRight;
			OneIndelEvent.WhetherReport = true;
         for (unsigned int GoodIndex = 1; GoodIndex < GoodNum; GoodIndex++) {
				//cout << GoodIndex << "\t" << GoodIndels[GoodIndex].BPLeft << "\t" << GoodIndels[GoodIndex].BPRight << "\t" << OneIndelEvent.BPLeft << "\t" << OneIndelEvent.BPRight << endl;
            if (GoodIndels[GoodIndex].BPLeft == OneIndelEvent.BPLeft
                && GoodIndels[GoodIndex].BPRight == OneIndelEvent.BPRight)
               OneIndelEvent.End = GoodIndex;
				else  {
					OneIndelEvent.RealStart = OneIndelEvent.BPLeft;
					OneIndelEvent.RealEnd = OneIndelEvent.BPRight;
					OneIndelEvent.Support = OneIndelEvent.End - OneIndelEvent.Start + 1;
					//GetRealStart4Deletion(CurrentChr, OneIndelEvent.RealStart, OneIndelEvent.RealEnd);
					IndelEvents.push_back(OneIndelEvent);
					OneIndelEvent.Start = GoodIndex;
					OneIndelEvent.End = GoodIndex;
					OneIndelEvent.BPLeft =  GoodIndels[GoodIndex].BPLeft;
					OneIndelEvent.BPRight =  GoodIndels[GoodIndex].BPRight;
            }
         }

         OneIndelEvent.RealStart = OneIndelEvent.BPLeft;
         OneIndelEvent.RealEnd = OneIndelEvent.BPRight;
         OneIndelEvent.Support = OneIndelEvent.End - OneIndelEvent.Start + 1;
			//cout << OneIndelEvent.Support << "\t" << OneIndelEvent.Start << "\t" << OneIndelEvent.End << endl;
         //GetRealStart4Deletion(CurrentChr, OneIndelEvent.RealStart, OneIndelEvent.RealEnd);
         IndelEvents.push_back(OneIndelEvent);
			//cout << "IndelEvent: " << IndelEvents.size() << endl;
			//string IndelType;
			unsigned int RealStart;
			unsigned int RealEnd;
			//bool WhetherDeletion = true;
			string IndelStr;
			unsigned int Max_Support;
			unsigned int Max_Support_Index;

			if (IndelEvents.size()) {
            for (unsigned EventIndex = 0; EventIndex < IndelEvents.size(); EventIndex++) {
					//cout << IndelEvents[EventIndex].Start << "\t" << IndelEvents[EventIndex].End << "\t" << IndelEvents[EventIndex].Support << endl;
					RealStart = IndelEvents[EventIndex].RealStart;
					RealEnd = IndelEvents[EventIndex].RealEnd;
					if (IndelEvents[EventIndex].Support < NumRead2ReportCutOff) continue;
					// report max one
					if (GoodIndels[IndelEvents[EventIndex].Start].IndelSize < BalanceCutoff) {
						OutputInversions(GoodIndels, CurrentChr,
											  IndelEvents[EventIndex].Start,
											  IndelEvents[EventIndex].End,
											  RealStart, RealEnd, InvOutf);
						NumberOfInvInstances++;
						Count_INV_NT_output++;
					}
					else if (ReportEvent(GoodIndels, IndelEvents[EventIndex].Start, IndelEvents[EventIndex].End)) {
						OutputInversions(GoodIndels, CurrentChr,
											  IndelEvents[EventIndex].Start,
											  IndelEvents[EventIndex].End,
											  RealStart, RealEnd, InvOutf);
						NumberOfInvInstances++;
						Count_INV_NT_output++;
					}
            }
         }
      }   // if (!Deletions[Box_index].empty())
   } // for (unsigned Box_index = 0; Box_index < NumBoxes; Box_index++)
   cout << "Inversions with non-template sequence (INV_NT): " << Count_INV_NT_output << endl << endl;
}

void SortOutputDI(const unsigned & NumBoxes, const string & CurrentChr, vector <SPLIT_READ> & Reads, vector <unsigned> DI[], ofstream & DIOutf) {
   cout << "Sorting and outputing deletions with non-template sequences ..." << endl;
   unsigned int DINum;
   short CompareResult;
   SPLIT_READ Temp4Exchange;
	/*
   unsigned int C_S = 0;
   unsigned int C_E = 0;
   unsigned int C_BP_Left;// = GoodSIs[0].BPLeft;
   unsigned int C_BP_Right;// = GoodSIs[0].BPRight;
   unsigned int C_Indelsize;
	 */
   unsigned int GoodNum;
	vector <SPLIT_READ> InputIndels;
   vector <SPLIT_READ> GoodIndels;
   vector <Indel4output> IndelEvents;

   for (unsigned Box_index = 0; Box_index < NumBoxes; Box_index++) {
      //cout << "Box_index: "   << Box_index << endl;
      if (!DI[Box_index].empty()) {
         DINum = DI[Box_index].size();
			InputIndels.clear();
			for (int i = 0; i < DINum; i++) InputIndels.push_back(Reads[DI[Box_index][i]]);
         for (unsigned int First = 0; First < DINum - 1; First++) {
            //if (InputIndels[First].OK)
				{
               for (unsigned int Second = First + 1; Second < DINum; Second++) {
                  //if (InputIndels[Second].OK)
						{
							if (InputIndels[First].BPLeft < InputIndels[Second].BPLeft) continue;
                     else if (InputIndels[First].BPLeft > InputIndels[Second].BPLeft) {
								CompareResult = 1;
							}
							else if (InputIndels[First].BPLeft == InputIndels[Second].BPLeft) {
								if (InputIndels[First].BPRight < InputIndels[Second].BPRight) continue;
								else if (InputIndels[First].BPRight > InputIndels[Second].BPRight) {
									CompareResult = 1;
								}
								else {
									if (InputIndels[First].NT_size < InputIndels[Second].NT_size) continue;
									else if (InputIndels[First].NT_size > InputIndels[Second].NT_size) CompareResult = 1;
                                    else if (CompareTwoString(InputIndels[First].NT_str, InputIndels[Second].NT_str)) { // NT_size ==
                                        CompareResult = 1;
                                    }
									else CompareResult = 2;
									//else {
									//	if (InputIndels[First].MatchedRelPos == InputIndels[Second].MatchedRelPos) {
									//		if (InputIndels[First].UnmatchedSeq == InputIndels[Second].UnmatchedSeq) {
									//			InputIndels[Second].OK = false;
									//		}
									//
									//	}
									//}
							   }
                     }
							if (CompareResult == 1) {
								Temp4Exchange = InputIndels[First];
  	                     InputIndels[First] = InputIndels[Second];
  	                     InputIndels[Second] = Temp4Exchange;
							}
							else if (CompareResult == 2) {
								Temp4Exchange = InputIndels[First + 1];
  	                     InputIndels[First + 1] = InputIndels[Second];
  	                     InputIndels[Second] = Temp4Exchange;
							}
						}
               }
            }
         }
         GoodIndels.clear();
	      IndelEvents.clear();

         for (unsigned int First = 0; First < DINum; First++) {
            //if (InputIndels[First].OK)
					GoodIndels.push_back(InputIndels[First]);
         }
         GoodNum = GoodIndels.size();
			//cout << Box_index << " " << GoodNum << endl;
			if (GoodNum == 0) continue;
         //cout << "GoodNum: " << GoodNum << endl;   string InsertedStr;   string NT_str;  short NT_size;
         Indel4output OneIndelEvent;
         OneIndelEvent.Start = 0;
	 OneIndelEvent.End = 0;
	 OneIndelEvent.IndelSize = GoodIndels[0].IndelSize;
         OneIndelEvent.NT_size = GoodIndels[0].NT_size;
	 OneIndelEvent.BPLeft =  GoodIndels[0].BPLeft;
	 OneIndelEvent.BPRight =  GoodIndels[0].BPRight;
     OneIndelEvent.IndelStr = GoodIndels[0].NT_str;
	 OneIndelEvent.WhetherReport = true;
         for (unsigned int GoodIndex = 1; GoodIndex < GoodNum; GoodIndex++) {
            if (GoodIndels[GoodIndex].BPLeft == OneIndelEvent.BPLeft
                && GoodIndels[GoodIndex].IndelSize == OneIndelEvent.IndelSize
                && GoodIndels[GoodIndex].NT_size == OneIndelEvent.NT_size
                && OneIndelEvent.IndelStr == GoodIndels[GoodIndex].NT_str
                )
               OneIndelEvent.End = GoodIndex;
	    else  {
	           IndelEvents.push_back(OneIndelEvent);
                   OneIndelEvent.Start = GoodIndex;
		   OneIndelEvent.End = GoodIndex;
		   OneIndelEvent.BPLeft =  GoodIndels[GoodIndex].BPLeft;
		   OneIndelEvent.IndelSize =  GoodIndels[GoodIndex].IndelSize;
                   OneIndelEvent.NT_size = GoodIndels[GoodIndex].NT_size;
                   OneIndelEvent.IndelStr = GoodIndels[GoodIndex].NT_str;
            }
         }

	     //if (OneIndelEvent.End - OneIndelEvent.Start + 1 >= NumRead2ReportCutOff)
	     IndelEvents.push_back(OneIndelEvent);
             unsigned int RealStart;
             unsigned int RealEnd;
             for (unsigned EventIndex = 0; EventIndex < IndelEvents.size(); EventIndex++) {
                  if (IndelEvents[EventIndex].End - IndelEvents[EventIndex].Start + 1 >= NumRead2ReportCutOff)
		  {
                     RealStart = GoodIndels[IndelEvents[EventIndex].Start].BPLeft;
                     RealEnd = GoodIndels[IndelEvents[EventIndex].Start].BPRight;
                     //if (IndelEvents[EventIndex].IndelSize < 100) {}
		     //if (IndelSize < BalanceCutoff) {
		     //        OutputDI(GoodIndels, CurrentChr,
		     //                 IndelEvents[EventIndex].Start,
		     //       	        IndelEvents[EventIndex].End,
		     // 		RealStart, RealEnd, DIOutf);
		     //        NumberOfDIInstances++;
	             //}
			   if (GoodIndels[IndelEvents[EventIndex].Start].IndelSize < BalanceCutoff) {
					OutputDI(GoodIndels, CurrentChr,
								IndelEvents[EventIndex].Start,
								IndelEvents[EventIndex].End,
								RealStart, RealEnd, DIOutf);
					NumberOfDIInstances++;
				}
		      else if (ReportEvent(GoodIndels, IndelEvents[EventIndex].Start, IndelEvents[EventIndex].End)) {
		              OutputDI(GoodIndels, CurrentChr,
		                       IndelEvents[EventIndex].Start,
		                       IndelEvents[EventIndex].End,
		                       RealStart, RealEnd, DIOutf);
		             NumberOfDIInstances++;
		     }
		  }
             }
      }   // if (!Deletions[Box_index].empty())
   } // for (unsigned Box_index = 0; Box_index < NumBoxes; Box_index++)
   cout << "deletions with non-template sequences: " << NumberOfDIInstances << endl << endl;
}


bool NotInVector(const string & OneTag, const vector <string> & VectorTag) {
  for (unsigned i = 0; i < VectorTag.size(); i++) {
     if (OneTag == VectorTag[i]) return false;
  }
  return true;
}


bool ReportEvent(const vector <SPLIT_READ> & Deletions, const unsigned int & S, const unsigned int & E) {
   short ReadLength = Deletions[S].ReadLength;
   short Min_Length = (short)((ReadLength * Min_Filter_Ratio) + 0.5) - 1;
   short Max_Length = (short)(ReadLength * (1 - Min_Filter_Ratio) - 0.5) - 1;
   bool LeftMin = false;
   bool LeftMax = false;
	bool RightMin = false;
   bool RightMax = false;
   for (unsigned i = S; i <= E; i++) {
		ReadLength = Deletions[i].ReadLength;
		Min_Length = (short)((ReadLength * Min_Filter_Ratio) + 0.5) - 1;
		Max_Length = (short)(ReadLength * (1 - Min_Filter_Ratio) - 0.5) - 1;
      if (Deletions[i].BP <= Min_Length) {
         LeftMin = true;
	      //break;
      }
		if (Deletions[i].ReadLength - Deletions[i].BP - Deletions[i].NT_size <= Min_Length) {
         RightMin = true;
	      //break;
      }
		if (Deletions[i].BP >= Max_Length) {
         LeftMax = true;
	      //break;
      }
		if (Deletions[i].ReadLength - Deletions[i].BP - Deletions[i].NT_size >= Max_Length) {
         RightMax = true;
	      //break;
      }
   }

   if (LeftMin && LeftMax && RightMin && RightMax) return true;
   else return false;
}

void GetRealStart4Deletion(const string & TheInput, unsigned int & RealStart, unsigned int & RealEnd) {
   unsigned int PosIndex = RealStart + SpacerBeforeAfter;
   unsigned int Start = PosIndex + 1;
   unsigned int End = RealEnd + SpacerBeforeAfter - 1;
   while (TheInput[PosIndex] == TheInput[End]) {
      --PosIndex;
      --End;
   }
   RealStart = PosIndex - SpacerBeforeAfter;
   PosIndex = RealEnd + SpacerBeforeAfter;
   while (TheInput[PosIndex] == TheInput[Start]) {
      ++PosIndex;
      ++Start;
   }
   RealEnd = PosIndex - SpacerBeforeAfter;
}

void GetRealStart4Insertion(const string & TheInput,
                            string & InsertedStr,
                            unsigned int &RealStart,
                            unsigned int &RealEnd)
{
    unsigned int IndelSize = InsertedStr.size();
    unsigned int PosIndex = RealStart + SpacerBeforeAfter;
    unsigned int original_RealStart = RealStart;

    for (int i = IndelSize - 1; i >= 0; i--) {
        if (TheInput[PosIndex] == InsertedStr[i]) {
            PosIndex--;
        }
        else {
            break;
        }
    }
    if (PosIndex == RealStart + SpacerBeforeAfter - IndelSize) {
        while (TheInput[PosIndex] == TheInput[PosIndex + IndelSize]) {
            PosIndex--;
        }
    }
    RealStart = PosIndex - SpacerBeforeAfter;
    PosIndex = RealEnd + SpacerBeforeAfter;
    for (unsigned int i = 0; i < IndelSize; i++) {
        if (TheInput[PosIndex] == InsertedStr[i]) {
            PosIndex++;
        }
        else {
            break;
        }
    }
    if (PosIndex == RealEnd + SpacerBeforeAfter + IndelSize) {
        while (TheInput[PosIndex] == TheInput[PosIndex - IndelSize]) {
            PosIndex++;
        }
    }
    RealEnd = PosIndex - SpacerBeforeAfter;
    unsigned DIFF = RealStart - original_RealStart;
    InsertedStr = InsertedStr.substr(0, IndelSize - DIFF) + InsertedStr.substr(IndelSize, DIFF);
}

vector <Region> Merge(const vector <Region> & AllRegions) {
	return AllRegions;
}

void GetCloseEnd(const string & CurrentChr, SPLIT_READ & Temp_One_Read) {

	Temp_One_Read.ReadLength = Temp_One_Read.UnmatchedSeq.size();
	Temp_One_Read.ReadLengthMinus = Temp_One_Read.ReadLength - 1;
	char LeftChar, RightChar;
	string CurrentReadSeq;
	vector <unsigned int> PD[Temp_One_Read.TOTAL_SNP_ERROR_CHECKED];
	for (int CheckIndex = 0; CheckIndex < Temp_One_Read.TOTAL_SNP_ERROR_CHECKED; CheckIndex++) {
		PD[CheckIndex].reserve(3 * Temp_One_Read.InsertSize);
	}
	vector <UniquePoint> UP;
	//char Direction;
	int Start, End;
	short BP_Start;// = MinClose;
	short BP_End;// = ReadLength - MinClose;

	//for (int i = 0; i < TOTAL_SNP_ERROR_CHECKED; i++) {
	//	PD[i].clear();
	//}
	//UP.clear();
	Temp_One_Read.UP_Close.clear();
	//MinClose = short(log((double)(Temp_One_Read.InsertSize * 3))/log(4.0) + 0.8) + 3;
	BP_Start = Temp_One_Read.MinClose;
	BP_End = Temp_One_Read.ReadLengthMinus;
	//Temp_One_Read.OK = true;
	if (Temp_One_Read.TOTAL_SNP_ERROR_CHECKED_Minus) {
		if (Temp_One_Read.MatchedD == Plus) {
			CurrentReadSeq = ReverseComplement(Temp_One_Read.UnmatchedSeq);
			Start = Temp_One_Read.MatchedRelPos + SpacerBeforeAfter;
			End = Start + 3 * Temp_One_Read.InsertSize;
			LeftChar = CurrentReadSeq[0];
			if (LeftChar != 'N') {
				for (int pos = Start; pos < End; pos++) {
					if (CurrentChr[pos] == LeftChar) {
						PD[0].push_back(pos);
					}
					else PD[1].push_back(pos);
				}
			}
			else { //Match2N[(short)'A'] = 'N';
				for (int pos = Start; pos < End; pos++) {
					if (Match2N[(short)CurrentChr[pos]] == 'N') {
						PD[0].push_back(pos);
					}
					//else Left_PD[1].push_back(pos);
				}
			}
			CheckLeft_Close(Temp_One_Read, CurrentChr, CurrentReadSeq, PD, BP_Start, BP_End, FirstBase, UP);
			//Direction = Minus;
			for (unsigned LeftUP_index = 0; LeftUP_index < UP.size(); LeftUP_index++) {
				if (CheckMismatches(CurrentChr, Temp_One_Read.UnmatchedSeq, UP[LeftUP_index])) {
					Temp_One_Read.Used = false;
					Temp_One_Read.UP_Close.push_back(UP[LeftUP_index]);
				}
			}
			UP.clear();
		}
		else if (Temp_One_Read.MatchedD == Minus) {
			CurrentReadSeq = Temp_One_Read.UnmatchedSeq;
			End = Temp_One_Read.MatchedRelPos + SpacerBeforeAfter;
			Start = End - 3 * Temp_One_Read.InsertSize;
			RightChar = CurrentReadSeq[Temp_One_Read.ReadLengthMinus];
			if (RightChar != 'N') {
				for (int pos = Start; pos < End; pos++) {
					if (CurrentChr[pos] == RightChar) {
						PD[0].push_back(pos);
					}
					else PD[1].push_back(pos);
				}
			}
			else { //Match2N[(short)'A'] = 'N';
				for (int pos = Start; pos < End; pos++) {
					if (Match2N[(short)CurrentChr[pos]] == 'N') {
						PD[0].push_back(pos);
					}
					//else Left_PD[1].push_back(pos);
				}
			}
			//cout << "1\t" << PD[0].size() << "\t" << PD[1].size() << endl;
			CheckRight_Close(Temp_One_Read, CurrentChr, CurrentReadSeq, PD, BP_Start, BP_End, FirstBase, UP);
			//cout << UP.size() << endl;
			//Direction = '+';
			for (unsigned RightUP_index = 0; RightUP_index < UP.size(); RightUP_index++) {
				if (CheckMismatches(CurrentChr, Temp_One_Read.UnmatchedSeq, UP[RightUP_index])) {
					Temp_One_Read.Used = false;
					Temp_One_Read.UP_Close.push_back(UP[RightUP_index]);
				}
			}
			UP.clear();
		}
	}
	else { // TOTAL_SNP_ERROR_CHECKED_Minus
		if (Temp_One_Read.MatchedD == Plus) {
			CurrentReadSeq = ReverseComplement(Temp_One_Read.UnmatchedSeq);
			Start = Temp_One_Read.MatchedRelPos + SpacerBeforeAfter;
			End = Start + 3 * Temp_One_Read.InsertSize;
			LeftChar = CurrentReadSeq[0];
			if (LeftChar != 'N') {
				for (int pos = Start; pos < End; pos++) {
					if (CurrentChr[pos] == LeftChar) {
						PD[0].push_back(pos);
					}
					//else PD[1].push_back(pos);
				}
			}
			else { //Match2N[(short)'A'] = 'N';
				for (int pos = Start; pos < End; pos++) {
					if (Match2N[(short)CurrentChr[pos]] == 'N') {
						PD[0].push_back(pos);
					}
					//else Left_PD[1].push_back(pos);
				}
			}
			CheckLeft_Close(Temp_One_Read, CurrentChr, CurrentReadSeq, PD, BP_Start, BP_End, FirstBase, UP);
			//Direction = Minus;
			for (unsigned LeftUP_index = 0; LeftUP_index < UP.size(); LeftUP_index++) {
				if (CheckMismatches(CurrentChr, Temp_One_Read.UnmatchedSeq, UP[LeftUP_index])) {
					Temp_One_Read.Used = false;
					Temp_One_Read.UP_Close.push_back(UP[LeftUP_index]);
				}
			}
			UP.clear();
		}
		else if (Temp_One_Read.MatchedD == Minus) {
			CurrentReadSeq = Temp_One_Read.UnmatchedSeq;
			End = Temp_One_Read.MatchedRelPos + SpacerBeforeAfter;
			Start = End - 3 * Temp_One_Read.InsertSize;
			RightChar = CurrentReadSeq[Temp_One_Read.ReadLengthMinus];
			if (RightChar != 'N') {
				for (int pos = Start; pos < End; pos++) {
					if (CurrentChr[pos] == RightChar) {
						PD[0].push_back(pos);
					}
					//else PD[1].push_back(pos);
				}
			}
			else { //Match2N[(short)'A'] = 'N';
				for (int pos = Start; pos < End; pos++) {
					if (Match2N[(short)CurrentChr[pos]] == 'N') {
						PD[0].push_back(pos);
					}
					//else Left_PD[1].push_back(pos);
				}
			}
			//cout << "2" << PD[0].size() << "\t" << PD[1].size() << endl;
			CheckRight_Close(Temp_One_Read, CurrentChr, CurrentReadSeq, PD, BP_Start, BP_End, FirstBase, UP);
			//cout << UP.size() << endl;
			//Direction = '+';
			for (unsigned RightUP_index = 0; RightUP_index < UP.size(); RightUP_index++) {
				if (CheckMismatches(CurrentChr, Temp_One_Read.UnmatchedSeq, UP[RightUP_index])) {
					Temp_One_Read.Used = false;
					Temp_One_Read.UP_Close.push_back(UP[RightUP_index]);
				}
			}
			UP.clear();
		}
	}
}

void GetFarEnd_SingleStrandDownStreamInsertions(const string & CurrentChr, SPLIT_READ & Temp_One_Read, const short & RangeIndex) {
	Temp_One_Read.ReadLength = Temp_One_Read.UnmatchedSeq.size();
	Temp_One_Read.ReadLengthMinus = Temp_One_Read.ReadLength - 1;
	char LeftChar, RightChar;
	string CurrentReadSeq;
	vector <unsigned int> PD[Temp_One_Read.TOTAL_SNP_ERROR_CHECKED];
	for (int CheckIndex = 0; CheckIndex < Temp_One_Read.TOTAL_SNP_ERROR_CHECKED; CheckIndex++) {
		PD[CheckIndex].reserve(Temp_One_Read.InsertSize * 2 + DSizeArray[RangeIndex]);
		//PD_Minus[CheckIndex].reserve(End - Start + 1);
	}
	vector <UniquePoint> UP;
	//char Direction;
	int Start, End;
	short BP_Start;// = MinClose;
	short BP_End;// = ReadLength - MinClose;

	//for (int i = 0; i < TOTAL_SNP_ERROR_CHECKED; i++) {
	//	PD[i].clear();
	//}
	//UP.clear();
	//Temp_One_Read.UP_Close.clear();
	//MinClose = short(log((double)(Temp_One_Read.InsertSize * 2 + Range))/log(4.0) + 0.8) + 3;
	BP_Start = Temp_One_Read.MinClose + RangeIndex;
	BP_End = Temp_One_Read.ReadLengthMinus;
	//Temp_One_Read.OK = true;
	if (Temp_One_Read.MatchedD == Minus) {
		//CurrentReadSeq = ReverseComplement(Temp_One_Read.UnmatchedSeq);
		CurrentReadSeq = Temp_One_Read.UnmatchedSeq;

		End = Temp_One_Read.UP_Close[0].AbsLoc + Temp_One_Read.UP_Close[0].LengthStr;
		if (End > SpacerBeforeAfter + Temp_One_Read.InsertSize * 2 + DSizeArray[RangeIndex])
		   Start = End - DSizeArray[RangeIndex] - Temp_One_Read.InsertSize * 2;
		else Start = SpacerBeforeAfter;


		if (End > CurrentChr.size() - SpacerBeforeAfter) End = CurrentChr.size() - SpacerBeforeAfter;
		LeftChar = CurrentReadSeq[0];
		if (Temp_One_Read.TOTAL_SNP_ERROR_CHECKED_Minus) {
			if (LeftChar != 'N') {
				for (int pos = Start; pos < End; pos++) {
					if (CurrentChr[pos] == LeftChar) {
						PD[0].push_back(pos);
					}
					else PD[1].push_back(pos);
				}
			}
			else { //Match2N[(short)'A'] = 'N';
				for (int pos = Start; pos < End; pos++) {
					if (Match2N[(short)CurrentChr[pos]] == 'N') {
						PD[0].push_back(pos);
					}
					//else Left_PD[1].push_back(pos);
				}
			}
		}
      else { // TOTAL_SNP_ERROR_CHECKED_Minus
			if (LeftChar != 'N') {
				for (int pos = Start; pos < End; pos++) {
					if (CurrentChr[pos] == LeftChar) {
						PD[0].push_back(pos);
					}
					//else PD[1].push_back(pos);
				}
			}
			else { //Match2N[(short)'A'] = 'N';
				for (int pos = Start; pos < End; pos++) {
					if (Match2N[(short)CurrentChr[pos]] == 'N') {
						PD[0].push_back(pos);
					}
					//else Left_PD[1].push_back(pos);
				}
			}
		}
		CheckLeft_Far(Temp_One_Read, CurrentChr, CurrentReadSeq, PD, BP_Start, BP_End, FirstBase, UP);
		//Direction = Minus;
		for (unsigned LeftUP_index = 0; LeftUP_index < UP.size(); LeftUP_index++) {
			if (CheckMismatches(CurrentChr, Temp_One_Read.UnmatchedSeq, UP[LeftUP_index])) {
				Temp_One_Read.Used = false;
				Temp_One_Read.UP_Far.push_back(UP[LeftUP_index]);
			}
		}
		UP.clear();
	}
	else if (Temp_One_Read.MatchedD == Plus) {
		CurrentReadSeq = ReverseComplement(Temp_One_Read.UnmatchedSeq);

		Start = Temp_One_Read.UP_Close[0].AbsLoc - Temp_One_Read.UP_Close[0].LengthStr;
		//Start = Temp_One_Read.MatchedRelPos + SpacerBeforeAfter;
		End = Start + DSizeArray[RangeIndex] + Temp_One_Read.InsertSize * 2;
		if (End > CurrentChr.size() - SpacerBeforeAfter) End = CurrentChr.size() - SpacerBeforeAfter;

		RightChar = CurrentReadSeq[Temp_One_Read.ReadLengthMinus];
		if (Temp_One_Read.TOTAL_SNP_ERROR_CHECKED_Minus) {
			if (RightChar != 'N') {
				for (int pos = Start; pos < End; pos++) {
					if (CurrentChr[pos] == RightChar) {
						PD[0].push_back(pos);
					}
					else PD[1].push_back(pos);
				}
			}
			else { //Match2N[(short)'A'] = 'N';
				for (int pos = Start; pos < End; pos++) {
					if (Match2N[(short)CurrentChr[pos]] == 'N') {
						PD[0].push_back(pos);
					}
					//else Left_PD[1].push_back(pos);
				}
			}
		}
		else {
			if (RightChar != 'N') {
				for (int pos = Start; pos < End; pos++) {
					if (CurrentChr[pos] == RightChar) {
						PD[0].push_back(pos);
					}
					//else PD[1].push_back(pos);
				}
			}
			else { //Match2N[(short)'A'] = 'N';
				for (int pos = Start; pos < End; pos++) {
					if (Match2N[(short)CurrentChr[pos]] == 'N') {
						PD[0].push_back(pos);
					}
					//else Left_PD[1].push_back(pos);
				}
			}
		}

		CheckRight_Far(Temp_One_Read, CurrentChr, CurrentReadSeq, PD, BP_Start, BP_End, FirstBase, UP);
		//Direction = '+';
		for (unsigned RightUP_index = 0; RightUP_index < UP.size(); RightUP_index++) {
			if (CheckMismatches(CurrentChr, Temp_One_Read.UnmatchedSeq, UP[RightUP_index])) {
				Temp_One_Read.Used = false;
				Temp_One_Read.UP_Far.push_back(UP[RightUP_index]);
			}
		}
		UP.clear();
	}
}

void GetFarEnd_SingleStrandDownStream(const string & CurrentChr, SPLIT_READ & Temp_One_Read, const short & RangeIndex) {
	Temp_One_Read.ReadLength = Temp_One_Read.UnmatchedSeq.size();
	Temp_One_Read.ReadLengthMinus = Temp_One_Read.ReadLength - 1;
	char LeftChar, RightChar;
	string CurrentReadSeq;
	vector <unsigned int> PD[Temp_One_Read.TOTAL_SNP_ERROR_CHECKED];
	for (int CheckIndex = 0; CheckIndex < Temp_One_Read.TOTAL_SNP_ERROR_CHECKED; CheckIndex++) {
		PD[CheckIndex].reserve(Temp_One_Read.InsertSize * 2 + DSizeArray[RangeIndex]);
		//PD_Minus[CheckIndex].reserve(End - Start + 1);
	}
	vector <UniquePoint> UP;
	//char Direction;
	int Start, End;
	short BP_Start;// = MinClose;
	short BP_End;// = ReadLength - MinClose;

	//for (int i = 0; i < TOTAL_SNP_ERROR_CHECKED; i++) {
	//	PD[i].clear();
	//}
	//UP.clear();
	//Temp_One_Read.UP_Close.clear();
	//MinClose = short(log((double)(Temp_One_Read.InsertSize * 2 + Range))/log(4.0) + 0.8) + 3;
	BP_Start = Temp_One_Read.MinClose + RangeIndex;
	BP_End = Temp_One_Read.ReadLengthMinus;
	//Temp_One_Read.OK = true;
	if (Temp_One_Read.MatchedD == Minus) {
		//CurrentReadSeq = ReverseComplement(Temp_One_Read.UnmatchedSeq);
		CurrentReadSeq = Temp_One_Read.UnmatchedSeq;

		End = Temp_One_Read.UP_Close[0].AbsLoc + Temp_One_Read.UP_Close[0].LengthStr - Temp_One_Read.ReadLength;
		if (End > SpacerBeforeAfter + Temp_One_Read.InsertSize * 2 + DSizeArray[RangeIndex])
		   Start = End - DSizeArray[RangeIndex] - Temp_One_Read.InsertSize * 2;
		else Start = SpacerBeforeAfter;


		if (End > CurrentChr.size() - SpacerBeforeAfter) End = CurrentChr.size() - SpacerBeforeAfter;
		LeftChar = CurrentReadSeq[0];
		if (Temp_One_Read.TOTAL_SNP_ERROR_CHECKED_Minus) {
			if (LeftChar != 'N') {
				for (int pos = Start; pos < End; pos++) {
					if (CurrentChr[pos] == LeftChar) {
						PD[0].push_back(pos);
					}
					else PD[1].push_back(pos);
				}
			}
			else { //Match2N[(short)'A'] = 'N';
				for (int pos = Start; pos < End; pos++) {
					if (Match2N[(short)CurrentChr[pos]] == 'N') {
						PD[0].push_back(pos);
					}
					//else Left_PD[1].push_back(pos);
				}
			}
		}
      else { // TOTAL_SNP_ERROR_CHECKED_Minus
			if (LeftChar != 'N') {
				for (int pos = Start; pos < End; pos++) {
					if (CurrentChr[pos] == LeftChar) {
						PD[0].push_back(pos);
					}
					//else PD[1].push_back(pos);
				}
			}
			else { //Match2N[(short)'A'] = 'N';
				for (int pos = Start; pos < End; pos++) {
					if (Match2N[(short)CurrentChr[pos]] == 'N') {
						PD[0].push_back(pos);
					}
					//else Left_PD[1].push_back(pos);
				}
			}
		}
		CheckLeft_Far(Temp_One_Read, CurrentChr, CurrentReadSeq, PD, BP_Start, BP_End, FirstBase, UP);
		//Direction = Minus;
		for (unsigned LeftUP_index = 0; LeftUP_index < UP.size(); LeftUP_index++) {
			if (CheckMismatches(CurrentChr, Temp_One_Read.UnmatchedSeq, UP[LeftUP_index])) {
				Temp_One_Read.Used = false;
				Temp_One_Read.UP_Far.push_back(UP[LeftUP_index]);
			}
		}
		UP.clear();
	}
	else if (Temp_One_Read.MatchedD == Plus) {
		CurrentReadSeq = ReverseComplement(Temp_One_Read.UnmatchedSeq);

		Start = Temp_One_Read.UP_Close[0].AbsLoc - Temp_One_Read.UP_Close[0].LengthStr + Temp_One_Read.ReadLength;
		//Start = Temp_One_Read.MatchedRelPos + SpacerBeforeAfter;
		End = Start + DSizeArray[RangeIndex] + Temp_One_Read.InsertSize * 2;
		if (End > CurrentChr.size() - SpacerBeforeAfter) End = CurrentChr.size() - SpacerBeforeAfter;

		RightChar = CurrentReadSeq[Temp_One_Read.ReadLengthMinus];
		if (Temp_One_Read.TOTAL_SNP_ERROR_CHECKED_Minus) {
			if (RightChar != 'N') {
				for (int pos = Start; pos < End; pos++) {
					if (CurrentChr[pos] == RightChar) {
						PD[0].push_back(pos);
					}
					else PD[1].push_back(pos);
				}
			}
			else { //Match2N[(short)'A'] = 'N';
				for (int pos = Start; pos < End; pos++) {
					if (Match2N[(short)CurrentChr[pos]] == 'N') {
						PD[0].push_back(pos);
					}
					//else Left_PD[1].push_back(pos);
				}
			}
		}
		else { // TOTAL_SNP_ERROR_CHECKED_Minus
			if (RightChar != 'N') {
				for (int pos = Start; pos < End; pos++) {
					if (CurrentChr[pos] == RightChar) {
						PD[0].push_back(pos);
					}
					//else PD[1].push_back(pos);
				}
			}
			else { //Match2N[(short)'A'] = 'N';
				for (int pos = Start; pos < End; pos++) {
					if (Match2N[(short)CurrentChr[pos]] == 'N') {
						PD[0].push_back(pos);
					}
					//else Left_PD[1].push_back(pos);
				}
			}
		}

		CheckRight_Far(Temp_One_Read, CurrentChr, CurrentReadSeq, PD, BP_Start, BP_End, FirstBase, UP);
		//Direction = '+';
		for (unsigned RightUP_index = 0; RightUP_index < UP.size(); RightUP_index++) {
			if (CheckMismatches(CurrentChr, Temp_One_Read.UnmatchedSeq, UP[RightUP_index])) {
				Temp_One_Read.Used = false;
				Temp_One_Read.UP_Far.push_back(UP[RightUP_index]);
			}
		}
		UP.clear();
	}
}

void GetFarEnd_SingleStrandUpStream(const string & CurrentChr, SPLIT_READ & Temp_One_Read, const short & RangeIndex) {
	Temp_One_Read.ReadLength = Temp_One_Read.UnmatchedSeq.size();
	Temp_One_Read.ReadLengthMinus = Temp_One_Read.ReadLength - 1;
	char LeftChar, RightChar;
	string CurrentReadSeq;
	vector <unsigned int> PD[Temp_One_Read.TOTAL_SNP_ERROR_CHECKED];
	for (int CheckIndex = 0; CheckIndex < Temp_One_Read.TOTAL_SNP_ERROR_CHECKED; CheckIndex++) {
		PD[CheckIndex].reserve(Temp_One_Read.InsertSize * 2 + DSizeArray[RangeIndex]);
		//PD_Minus[CheckIndex].reserve(End - Start + 1);
	}
	vector <UniquePoint> UP;
	//char Direction;
	int Start, End;
	short BP_Start;// = MinClose;
	short BP_End;// = ReadLength - MinClose;

	//for (int i = 0; i < TOTAL_SNP_ERROR_CHECKED; i++) {
	//	PD[i].clear();
	//}
	//UP.clear();
	//Temp_One_Read.UP_Close.clear();
	//MinClose = short(log((double)(Temp_One_Read.InsertSize * 2 + DSizeArray[RangeIndex]))/log(4.0) + 0.8) + 3;
	BP_Start = Temp_One_Read.MinClose + RangeIndex;
	BP_End = Temp_One_Read.ReadLengthMinus;
	//Temp_One_Read.OK = true;
	if (Temp_One_Read.MatchedD == Minus) {
		//CurrentReadSeq = ReverseComplement(Temp_One_Read.UnmatchedSeq);
		CurrentReadSeq = Temp_One_Read.UnmatchedSeq;

		Start = Temp_One_Read.UP_Close[0].AbsLoc + Temp_One_Read.UP_Close[0].LengthStr;
		End = Start + DSizeArray[RangeIndex] + Temp_One_Read.InsertSize * 2;
		if (End > CurrentChr.size() - SpacerBeforeAfter) End = CurrentChr.size() - SpacerBeforeAfter;

		LeftChar = CurrentReadSeq[0];
		if (Temp_One_Read.TOTAL_SNP_ERROR_CHECKED_Minus) {
			if (LeftChar != 'N') {
				for (int pos = Start; pos < End; pos++) {
					if (CurrentChr[pos] == LeftChar) {
						PD[0].push_back(pos);
					}
					else PD[1].push_back(pos);
				}
			}
			else { //Match2N[(short)'A'] = 'N';
				for (int pos = Start; pos < End; pos++) {
					if (Match2N[(short)CurrentChr[pos]] == 'N') {
						PD[0].push_back(pos);
					}
					//else Left_PD[1].push_back(pos);
				}
			}
		}
      else { // TOTAL_SNP_ERROR_CHECKED_Minus
			if (LeftChar != 'N') {
				for (int pos = Start; pos < End; pos++) {
					if (CurrentChr[pos] == LeftChar) {
						PD[0].push_back(pos);
					}
					//else PD[1].push_back(pos);
				}
			}
			else { //Match2N[(short)'A'] = 'N';
				for (int pos = Start; pos < End; pos++) {
					if (Match2N[(short)CurrentChr[pos]] == 'N') {
						PD[0].push_back(pos);
					}
					//else Left_PD[1].push_back(pos);
				}
			}
		}
		CheckLeft_Far(Temp_One_Read, CurrentChr, CurrentReadSeq, PD, BP_Start, BP_End, FirstBase, UP);
		//Direction = Minus;
		for (unsigned LeftUP_index = 0; LeftUP_index < UP.size(); LeftUP_index++) {
			if (CheckMismatches(CurrentChr, Temp_One_Read.UnmatchedSeq, UP[LeftUP_index])) {
				Temp_One_Read.Used = false;
				Temp_One_Read.UP_Far.push_back(UP[LeftUP_index]);
			}
		}
		UP.clear();
	}
	else if (Temp_One_Read.MatchedD == Plus) {
		CurrentReadSeq = ReverseComplement(Temp_One_Read.UnmatchedSeq);

		End = Temp_One_Read.UP_Close[0].AbsLoc - Temp_One_Read.UP_Close[0].LengthStr;
		//Start = Temp_One_Read.MatchedRelPos + SpacerBeforeAfter;

		if (End > DSizeArray[RangeIndex] + Temp_One_Read.InsertSize * 2 + SpacerBeforeAfter)
			Start = End - DSizeArray[RangeIndex] - Temp_One_Read.InsertSize * 2;
		else Start = SpacerBeforeAfter;

		RightChar = CurrentReadSeq[Temp_One_Read.ReadLengthMinus];
		if (Temp_One_Read.TOTAL_SNP_ERROR_CHECKED_Minus) {
			if (RightChar != 'N') {
				for (int pos = Start; pos < End; pos++) {
					if (CurrentChr[pos] == RightChar) {
						PD[0].push_back(pos);
					}
					else PD[1].push_back(pos);
				}
			}
			else { //Match2N[(short)'A'] = 'N';
				for (int pos = Start; pos < End; pos++) {
					if (Match2N[(short)CurrentChr[pos]] == 'N') {
						PD[0].push_back(pos);
					}
					//else Left_PD[1].push_back(pos);
				}
			}
		}
      else {
			if (RightChar != 'N') {
				for (int pos = Start; pos < End; pos++) {
					if (CurrentChr[pos] == RightChar) {
						PD[0].push_back(pos);
					}
					//else PD[1].push_back(pos);
				}
			}
			else { //Match2N[(short)'A'] = 'N';
				for (int pos = Start; pos < End; pos++) {
					if (Match2N[(short)CurrentChr[pos]] == 'N') {
						PD[0].push_back(pos);
					}
					//else Left_PD[1].push_back(pos);
				}
			}
		}
		CheckRight_Far(Temp_One_Read, CurrentChr, CurrentReadSeq, PD, BP_Start, BP_End, FirstBase, UP);
		//Direction = '+';
		for (unsigned RightUP_index = 0; RightUP_index < UP.size(); RightUP_index++) {
			if (CheckMismatches(CurrentChr, Temp_One_Read.UnmatchedSeq, UP[RightUP_index])) {
				Temp_One_Read.Used = false;
				Temp_One_Read.UP_Far.push_back(UP[RightUP_index]);
			}
		}
		UP.clear();
	}
}

void GetFarEnd_OtherStrand(const string & CurrentChr, SPLIT_READ & Temp_One_Read, const short & RangeIndex) {
	//short ReadLength = Temp_One_Read.UnmatchedSeq.size();
	//short ReadLengthMinus = ReadLength - 1;
	//MinClose = short(log((double)(Temp_One_Read.InsertSize * 5 + Range * 2))/log(4.0) + 0.8) + 3;// + MAX_SNP_ERROR;
	int Start, End;
	short BP_Start = Temp_One_Read.MinClose + RangeIndex;
	short BP_End = Temp_One_Read.ReadLengthMinus;
	//char Direction;

	vector <UniquePoint> UP;
	vector <unsigned int> PD_Plus[Temp_One_Read.TOTAL_SNP_ERROR_CHECKED];
	vector <unsigned int> PD_Minus[Temp_One_Read.TOTAL_SNP_ERROR_CHECKED];

	if (Temp_One_Read.MatchedD == Plus) {
		Start = Temp_One_Read.MatchedRelPos + SpacerBeforeAfter - Temp_One_Read.ReadLength - 2 * Temp_One_Read.InsertSize - DSizeArray[RangeIndex] ;
		End = Temp_One_Read.MatchedRelPos + SpacerBeforeAfter - Temp_One_Read.ReadLength + 3 * Temp_One_Read.InsertSize + DSizeArray[RangeIndex];
	}
	else {
		Start = Temp_One_Read.MatchedRelPos + SpacerBeforeAfter - 3 * Temp_One_Read.InsertSize + Temp_One_Read.ReadLength - DSizeArray[RangeIndex];
		End = Temp_One_Read.MatchedRelPos + SpacerBeforeAfter + 2 * Temp_One_Read.InsertSize + Temp_One_Read.ReadLength + DSizeArray[RangeIndex];
	}

	for (int CheckIndex = 0; CheckIndex < Temp_One_Read.TOTAL_SNP_ERROR_CHECKED; CheckIndex++) {
		PD_Plus[CheckIndex].reserve(End - Start + 1);
		PD_Minus[CheckIndex].reserve(End - Start + 1);
	}

	char CurrentBase = Temp_One_Read.UnmatchedSeq[0];
	char CurrentBaseRC = Convert2RC4N[(short)CurrentBase];
	if (Temp_One_Read.TOTAL_SNP_ERROR_CHECKED_Minus) {
		if (Temp_One_Read.MatchedD == Plus) {
			if (CurrentBase != 'N') {
				for (int pos = Start; pos < End; pos++) {
					if (CurrentChr[pos] == CurrentBase) PD_Plus[0].push_back(pos);
					else PD_Plus[1].push_back(pos);
					//if (CurrentChr[pos] == CurrentBaseRC) PD_Minus[0].push_back(pos);
					//else PD_Minus[1].push_back(pos);
				}
			}
			else { //Match2N[(short)'A'] = 'N';
				for (int pos = Start; pos < End; pos++) {
					if (Match2N[(short)CurrentChr[pos]] == 'N') {
						PD_Plus[0].push_back(pos);
						//PD_Minus[0].push_back(pos);
					}
					else {
						PD_Plus[1].push_back(pos);
						//PD_Minus[1].push_back(pos);
					}
				}
			}
		}
		else { // -
			if (CurrentBase != 'N') {
				for (int pos = Start; pos < End; pos++) {
					//if (CurrentChr[pos] == CurrentBase) PD_Plus[0].push_back(pos);
					//else PD_Plus[1].push_back(pos);
					if (CurrentChr[pos] == CurrentBaseRC) PD_Minus[0].push_back(pos);
					else PD_Minus[1].push_back(pos);
				}
			}
			else { //Match2N[(short)'A'] = 'N';
				for (int pos = Start; pos < End; pos++) {
					if (Match2N[(short)CurrentChr[pos]] == 'N') {
						//PD_Plus[0].push_back(pos);
						PD_Minus[0].push_back(pos);
					}
					else {
						//PD_Plus[1].push_back(pos);
						PD_Minus[1].push_back(pos);
					}
				}
			}
		}
	}
   else { // TOTAL_SNP_ERROR_CHECKED_Minus
		if (Temp_One_Read.MatchedD == Plus) {
			if (CurrentBase != 'N') {
				for (int pos = Start; pos < End; pos++) {
					if (CurrentChr[pos] == CurrentBase) PD_Plus[0].push_back(pos);
					//else PD_Plus[1].push_back(pos);
					//if (CurrentChr[pos] == CurrentBaseRC) PD_Minus[0].push_back(pos);
					//else PD_Minus[1].push_back(pos);
				}
			}
			else { //Match2N[(short)'A'] = 'N';
				for (int pos = Start; pos < End; pos++) {
					if (Match2N[(short)CurrentChr[pos]] == 'N') {
						PD_Plus[0].push_back(pos);
						//PD_Minus[0].push_back(pos);
					}
					//else {
					//	PD_Plus[1].push_back(pos);
						//PD_Minus[1].push_back(pos);
					//}
				}
			}
		}
		else { // -
			if (CurrentBase != 'N') {
				for (int pos = Start; pos < End; pos++) {
					//if (CurrentChr[pos] == CurrentBase) PD_Plus[0].push_back(pos);
					//else PD_Plus[1].push_back(pos);
					if (CurrentChr[pos] == CurrentBaseRC) PD_Minus[0].push_back(pos);
					//else PD_Minus[1].push_back(pos);
				}
			}
			else { //Match2N[(short)'A'] = 'N';
				for (int pos = Start; pos < End; pos++) {
					if (Match2N[(short)CurrentChr[pos]] == 'N') {
						//PD_Plus[0].push_back(pos);
						PD_Minus[0].push_back(pos);
					}
					//else {
						//PD_Plus[1].push_back(pos);
					//	PD_Minus[1].push_back(pos);
					//}
				}
			}
		}
	}

	CheckBoth(Temp_One_Read, CurrentChr, Temp_One_Read.UnmatchedSeq, PD_Plus, PD_Minus, BP_Start, BP_End, FirstBase, UP);
	for (unsigned UP_index = 0; UP_index < UP.size(); UP_index++) {
		if (UP[UP_index].Direction == Plus) {
			//Direction = Minus;
			if (CheckMismatches(CurrentChr, Temp_One_Read.UnmatchedSeq, UP[UP_index]))
				Temp_One_Read.UP_Far.push_back(UP[UP_index]);
		}
		else {
			//Direction = '+';
			if (CheckMismatches(CurrentChr, Temp_One_Read.UnmatchedSeq, UP[UP_index]))
				Temp_One_Read.UP_Far.push_back(UP[UP_index]);
		}
	}
	UP.clear();
	if (Temp_One_Read.UP_Far.size()) {
		//Before = Temp_One_Read.UP_Far.size();
		//cout << "1: " << Temp_One_Read.UP_Far.size() << "\tt\t";
		CleanUniquePoints(Temp_One_Read.UP_Far);
		//After = Temp_One_Read.UP_Far.size();
		//if (Before != After)
		//   cout << Before << "\t" << After << endl;
	}
	return;
}

void GetFarEnd_BothStrands(const string & CurrentChr, SPLIT_READ & Temp_One_Read, const short & RangeIndex) {
	//short ReadLength = Temp_One_Read.UnmatchedSeq.size();
	//short ReadLengthMinus = ReadLength - 1;
	//MinClose = short(log((double)(Temp_One_Read.InsertSize * 5 + DSizeArray[RangeIndex] * 2))/log(4.0) + 0.8) + 3;// + MAX_SNP_ERROR;
	int Start, End;
	short BP_Start = Temp_One_Read.MinClose + RangeIndex;
	short BP_End = Temp_One_Read.ReadLengthMinus;
	//char Direction;

	vector <UniquePoint> UP;
	vector <unsigned int> PD_Plus[Temp_One_Read.TOTAL_SNP_ERROR_CHECKED];
	vector <unsigned int> PD_Minus[Temp_One_Read.TOTAL_SNP_ERROR_CHECKED];



	if (Temp_One_Read.MatchedD == Plus) {
		Start = Temp_One_Read.MatchedRelPos + SpacerBeforeAfter - Temp_One_Read.ReadLength - 2 * Temp_One_Read.InsertSize - DSizeArray[RangeIndex] ;
		End = Temp_One_Read.MatchedRelPos + SpacerBeforeAfter - Temp_One_Read.ReadLength + 3 * Temp_One_Read.InsertSize + DSizeArray[RangeIndex];
	}
	else {
		Start = Temp_One_Read.MatchedRelPos + SpacerBeforeAfter - 3 * Temp_One_Read.InsertSize + Temp_One_Read.ReadLength - DSizeArray[RangeIndex];
		End = Temp_One_Read.MatchedRelPos + SpacerBeforeAfter + 2 * Temp_One_Read.InsertSize + Temp_One_Read.ReadLength + DSizeArray[RangeIndex];
	}

	for (int CheckIndex = 0; CheckIndex < Temp_One_Read.TOTAL_SNP_ERROR_CHECKED; CheckIndex++) {
		PD_Plus[CheckIndex].reserve(End - Start + 1);
		PD_Minus[CheckIndex].reserve(End - Start + 1);
	}

	char CurrentBase = Temp_One_Read.UnmatchedSeq[0];
	char CurrentBaseRC = Convert2RC4N[(short)CurrentBase];
	if (Temp_One_Read.TOTAL_SNP_ERROR_CHECKED_Minus) {
		if (CurrentBase != 'N') {
			for (int pos = Start; pos < End; pos++) {
				if (CurrentChr[pos] == CurrentBase) PD_Plus[0].push_back(pos);
				else PD_Plus[1].push_back(pos);
				if (CurrentChr[pos] == CurrentBaseRC) PD_Minus[0].push_back(pos);
				else PD_Minus[1].push_back(pos);
			}
		}
		else { //Match2N[(short)'A'] = 'N';
			for (int pos = Start; pos < End; pos++) {
				if (Match2N[(short)CurrentChr[pos]] == 'N') {
					PD_Plus[0].push_back(pos);
					PD_Minus[0].push_back(pos);
				}
				else {
					PD_Plus[1].push_back(pos);
					PD_Minus[1].push_back(pos);
				}
			}
		}
	}
	else { // TOTAL_SNP_ERROR_CHECKED_Minus
		if (CurrentBase != 'N') {
			for (int pos = Start; pos < End; pos++) {
				if (CurrentChr[pos] == CurrentBase) PD_Plus[0].push_back(pos);
				//else PD_Plus[1].push_back(pos);
				if (CurrentChr[pos] == CurrentBaseRC) PD_Minus[0].push_back(pos);
				//else PD_Minus[1].push_back(pos);
			}
		}
		else { //Match2N[(short)'A'] = 'N';
			for (int pos = Start; pos < End; pos++) {
				if (Match2N[(short)CurrentChr[pos]] == 'N') {
					PD_Plus[0].push_back(pos);
					PD_Minus[0].push_back(pos);
				}
				//else {
				//	PD_Plus[1].push_back(pos);
				//	PD_Minus[1].push_back(pos);
				//}
			}
		}
	}

	CheckBoth(Temp_One_Read, CurrentChr, Temp_One_Read.UnmatchedSeq, PD_Plus, PD_Minus, BP_Start, BP_End, FirstBase, UP);
	for (unsigned UP_index = 0; UP_index < UP.size(); UP_index++) {
		if (UP[UP_index].Direction == Plus) {
			//Direction = Minus;
			if (CheckMismatches(CurrentChr, Temp_One_Read.UnmatchedSeq, UP[UP_index]))
				Temp_One_Read.UP_Far.push_back(UP[UP_index]);
		}
		else {
			//Direction = '+';
			if (CheckMismatches(CurrentChr, Temp_One_Read.UnmatchedSeq, UP[UP_index]))
				Temp_One_Read.UP_Far.push_back(UP[UP_index]);
		}
	}
	UP.clear();
	if (Temp_One_Read.UP_Far.size()) {
		//Before = Temp_One_Read.UP_Far.size();
		//cout << "1: " << Temp_One_Read.UP_Far.size() << "\tt\t";
		CleanUniquePoints(Temp_One_Read.UP_Far);
		//After = Temp_One_Read.UP_Far.size();
		//if (Before != After)
		//   cout << Before << "\t" << After << endl;
	}
	return;
}

void GetFarEnd(const string & CurrentChr, SPLIT_READ & Temp_One_Read, const int & in_start, const int & in_end) {
	//short ReadLength = Temp_One_Read.UnmatchedSeq.size();
	//short ReadLengthMinus = ReadLength - 1;
	//MinClose = short(log((double)(in_end - in_start))/log(4.0) + 0.8) + 3 + MAX_SNP_ERROR;
	int Start, End;
	short BP_Start = Temp_One_Read.MinClose;
	short BP_End = Temp_One_Read.ReadLengthMinus;
	//char Direction;
	vector <UniquePoint> UP;
	vector <unsigned int> PD_Plus[Temp_One_Read.TOTAL_SNP_ERROR_CHECKED];
	vector <unsigned int> PD_Minus[Temp_One_Read.TOTAL_SNP_ERROR_CHECKED];

	for (int CheckIndex = 0; CheckIndex < Temp_One_Read.TOTAL_SNP_ERROR_CHECKED; CheckIndex++) {
		PD_Plus[CheckIndex].reserve(in_end - in_start + 1);
		PD_Minus[CheckIndex].reserve(in_end - in_start + 1);
	}
		Start = in_start;
		End = in_end;

	char CurrentBase = Temp_One_Read.UnmatchedSeq[0];
	char CurrentBaseRC = Convert2RC4N[(short)CurrentBase];
	if (Temp_One_Read.TOTAL_SNP_ERROR_CHECKED_Minus) {
		if (CurrentBase != 'N') {
			for (int pos = Start; pos < End; pos++) {
				//if (BreakDancerMask[pos] == BD_char) {
				if (CurrentChr[pos] == CurrentBase) PD_Plus[0].push_back(pos);
				else PD_Plus[1].push_back(pos);
				if (CurrentChr[pos] == CurrentBaseRC) PD_Minus[0].push_back(pos);
				else PD_Minus[1].push_back(pos);
				//}
			}
		}
		else { //Match2N[(short)'A'] = 'N';
			for (int pos = Start; pos < End; pos++) {
				//if (BreakDancerMask[pos] == BD_char) {
				if (Match2N[(short)CurrentChr[pos]] == 'N') {
					PD_Plus[0].push_back(pos);
					PD_Minus[0].push_back(pos);
				}
				else {
					PD_Plus[1].push_back(pos);
					PD_Minus[1].push_back(pos);
				}
				//}
			}
		}
	}
   else {
		if (CurrentBase != 'N') {
			for (int pos = Start; pos < End; pos++) {
				//if (BreakDancerMask[pos] == BD_char) {
				if (CurrentChr[pos] == CurrentBase) PD_Plus[0].push_back(pos);
				//else PD_Plus[1].push_back(pos);
				if (CurrentChr[pos] == CurrentBaseRC) PD_Minus[0].push_back(pos);
				//else PD_Minus[1].push_back(pos);
				//}
			}
		}
		else { //Match2N[(short)'A'] = 'N';
			for (int pos = Start; pos < End; pos++) {
				//if (BreakDancerMask[pos] == BD_char) {
				if (Match2N[(short)CurrentChr[pos]] == 'N') {
					PD_Plus[0].push_back(pos);
					PD_Minus[0].push_back(pos);
				}
				//else {
				//	PD_Plus[1].push_back(pos);
				//	PD_Minus[1].push_back(pos);
				//}
				//}
			}
		}
	}
	//cout << PD_Plus[0].size() << "\t" << PD_Plus[1].size() << endl;
	//cout << PD_Minus[0].size() << "\t" << PD_Minus[1].size() << endl;
	CheckBoth(Temp_One_Read, CurrentChr, Temp_One_Read.UnmatchedSeq, PD_Plus, PD_Minus, BP_Start, BP_End, FirstBase, UP);
	//cout << UP.size() << endl;
	for (unsigned UP_index = 0; UP_index < UP.size(); UP_index++) {
		if (UP[UP_index].Direction == Plus) {
			//Direction = Minus;
			if (CheckMismatches(CurrentChr, Temp_One_Read.UnmatchedSeq, UP[UP_index]))
				Temp_One_Read.UP_Far.push_back(UP[UP_index]);
		}
		else {
			//Direction = '+';
			if (CheckMismatches(CurrentChr, Temp_One_Read.UnmatchedSeq, UP[UP_index]))
				Temp_One_Read.UP_Far.push_back(UP[UP_index]);
		}
	}
	UP.clear();

	if (Temp_One_Read.UP_Far.size()) {
		//Before = Temp_One_Read.UP_Far.size();
		//cout << "1: " << Temp_One_Read.UP_Far.size() << "\tt\t";
		CleanUniquePoints(Temp_One_Read.UP_Far);
		//After = Temp_One_Read.UP_Far.size();
		//if (Before != After)
		//   cout << Before << "\t" << After << endl;
	}

	return;
}

void CheckBoth(const SPLIT_READ & OneRead,
					const string & TheInput,
               const string & CurrentReadSeq,
               const vector <unsigned int> PD_Plus[],
					const vector <unsigned int> PD_Minus[],
               const short & BP_Start,
               const short & BP_End,
               const short & CurrentLength,
               vector <UniquePoint> & UP) {
	int Sum;
   if (CurrentLength >= BP_Start && CurrentLength <= BP_End) {
      // put it to LeftUP if unique
		for (short i = 0; i <= OneRead.MAX_SNP_ERROR; i++) {
			if (PD_Plus[i].size() + PD_Minus[i].size() == 1 && CurrentLength >= BP_Start + i) {
				Sum = 0;
				if (ADDITIONAL_MISMATCH)
		   		for (short j = 1; j <= ADDITIONAL_MISMATCH; j++)
			   		Sum += PD_Plus[i + j].size() + PD_Minus[i + j].size();

				if (Sum == 0 && i <= (short)(Seq_Error_Rate * CurrentLength + 1)) {
					UniquePoint TempOne;
					TempOne.LengthStr = CurrentLength;
					if (PD_Plus[i].size() == 1) {
						TempOne.Direction = FORWARD;
						TempOne.Strand = SENSE;
						TempOne.AbsLoc = PD_Plus[i][0];
					}
					else if (PD_Minus[i].size() == 1) {
						TempOne.Direction = BACKWARD;
						TempOne.Strand = ANTISENSE;
						TempOne.AbsLoc = PD_Minus[i][0];
					}
					TempOne.Mismatches = i;
					UP.push_back(TempOne);
					break;
				}
			}
		}
   }
   if (CurrentLength < BP_End) {
      vector <unsigned int> PD_Plus_Output[OneRead.TOTAL_SNP_ERROR_CHECKED];
		vector <unsigned int> PD_Minus_Output[OneRead.TOTAL_SNP_ERROR_CHECKED];
		for (int CheckedIndex = 0; CheckedIndex < OneRead.TOTAL_SNP_ERROR_CHECKED; CheckedIndex++) {
			PD_Plus_Output[CheckedIndex].reserve(PD_Plus[CheckedIndex].size());
			PD_Minus_Output[CheckedIndex].reserve(PD_Minus[CheckedIndex].size());
		}
      const char CurrentChar = CurrentReadSeq[CurrentLength];
		const char CurrentCharRC = Convert2RC4N[(short)CurrentChar];
      //const int SizeOfCurrent = Left_PD.size();
      unsigned int pos;
		int SizeOfCurrent;
		//if (TOTAL_SNP_ERROR_CHECKED_Minus)
		{
			for (int i = 0; i < OneRead.TOTAL_SNP_ERROR_CHECKED_Minus; i++) {
				SizeOfCurrent = PD_Plus[i].size();
				if (CurrentChar == 'N') {
					for (int j = 0; j < SizeOfCurrent; j++) {
						pos =  PD_Plus[i][j] + 1;
						if (Match2N[(short)TheInput[pos]] == 'N')
							PD_Plus_Output[i].push_back(pos);
						else PD_Plus_Output[i + 1].push_back(pos);
					}
				}
				else {
					for (int j = 0; j < SizeOfCurrent; j++) {
						pos =  PD_Plus[i][j] + 1;
						if (TheInput[pos] == CurrentChar)
							PD_Plus_Output[i].push_back(pos);
						else PD_Plus_Output[i + 1].push_back(pos);
					}
				}
				SizeOfCurrent = PD_Minus[i].size();
				if (CurrentCharRC == 'N') {
					for (int j = 0; j < SizeOfCurrent; j++) {
						pos =  PD_Minus[i][j] - 1;
						if (Match2N[(short)TheInput[pos]] == 'N')
							PD_Minus_Output[i].push_back(pos);
						else PD_Minus_Output[i + 1].push_back(pos);
					}
				}
				else {
					for (int j = 0; j < SizeOfCurrent; j++) {
						pos =  PD_Minus[i][j] - 1;
						if (TheInput[pos] == CurrentCharRC)
							PD_Minus_Output[i].push_back(pos);
						else PD_Minus_Output[i + 1].push_back(pos);
					}
				}
			}

			SizeOfCurrent = PD_Plus[OneRead.TOTAL_SNP_ERROR_CHECKED_Minus].size();
			if (CurrentChar == 'N') {
				for (int j = 0; j < SizeOfCurrent; j++) {
					pos =  PD_Plus[OneRead.TOTAL_SNP_ERROR_CHECKED_Minus][j] + 1;
					if (Match2N[(short)TheInput[pos]] == 'N')
						PD_Plus_Output[OneRead.TOTAL_SNP_ERROR_CHECKED_Minus].push_back(pos);
				}
			}
			else {
				for (int j = 0; j < SizeOfCurrent; j++) {
					pos =  PD_Plus[OneRead.TOTAL_SNP_ERROR_CHECKED_Minus][j] + 1;
					if (TheInput[pos] == CurrentChar)
						PD_Plus_Output[OneRead.TOTAL_SNP_ERROR_CHECKED_Minus].push_back(pos);
				}
			}
			SizeOfCurrent = PD_Minus[OneRead.TOTAL_SNP_ERROR_CHECKED_Minus].size();
			if (CurrentCharRC == 'N') {
				for (int j = 0; j < SizeOfCurrent; j++) {
					pos =  PD_Minus[OneRead.TOTAL_SNP_ERROR_CHECKED_Minus][j] - 1;
					if (Match2N[(short)TheInput[pos]] == 'N')
						PD_Minus_Output[OneRead.TOTAL_SNP_ERROR_CHECKED_Minus].push_back(pos);
				}
			}
			else {
				for (int j = 0; j < SizeOfCurrent; j++) {
					pos =  PD_Minus[OneRead.TOTAL_SNP_ERROR_CHECKED_Minus][j] - 1;
					if (TheInput[pos] == CurrentCharRC)
						PD_Minus_Output[OneRead.TOTAL_SNP_ERROR_CHECKED_Minus].push_back(pos);
				}
			}


			Sum = 0;
			for (int i = 0; i <= OneRead.MAX_SNP_ERROR; i++) {
				Sum += PD_Plus_Output[i].size() + PD_Minus_Output[i].size();
			}
			if (Sum) {
				const short CurrentLengthOutput = CurrentLength + 1;
				CheckBoth(OneRead, TheInput, CurrentReadSeq, PD_Plus_Output, PD_Minus_Output,
							 BP_Start, BP_End,
							 CurrentLengthOutput, UP);
			}
			else return;
		}
		/*
      else { // TOTAL_SNP_ERROR_CHECKED_Minus
			SizeOfCurrent = PD_Plus[TOTAL_SNP_ERROR_CHECKED_Minus].size();
			if (CurrentChar == 'N') {
				for (int j = 0; j < SizeOfCurrent; j++) {
					pos =  PD_Plus[TOTAL_SNP_ERROR_CHECKED_Minus][j] + 1;
					if (Match2N[(short)TheInput[pos]] == 'N')
						PD_Plus_Output[TOTAL_SNP_ERROR_CHECKED_Minus].push_back(pos);
				}
			}
			else {
				for (int j = 0; j < SizeOfCurrent; j++) {
					pos =  PD_Plus[TOTAL_SNP_ERROR_CHECKED_Minus][j] + 1;
					if (TheInput[pos] == CurrentChar)
						PD_Plus_Output[TOTAL_SNP_ERROR_CHECKED_Minus].push_back(pos);
				}
			}
			SizeOfCurrent = PD_Minus[TOTAL_SNP_ERROR_CHECKED_Minus].size();
			if (CurrentCharRC == 'N') {
				for (int j = 0; j < SizeOfCurrent; j++) {
					pos =  PD_Minus[TOTAL_SNP_ERROR_CHECKED_Minus][j] - 1;
					if (Match2N[(short)TheInput[pos]] == 'N')
						PD_Minus_Output[TOTAL_SNP_ERROR_CHECKED_Minus].push_back(pos);
				}
			}
			else {
				for (int j = 0; j < SizeOfCurrent; j++) {
					pos =  PD_Minus[TOTAL_SNP_ERROR_CHECKED_Minus][j] - 1;
					if (TheInput[pos] == CurrentCharRC)
						PD_Minus_Output[TOTAL_SNP_ERROR_CHECKED_Minus].push_back(pos);
				}
			}

			if (PD_Plus_Output[TOTAL_SNP_ERROR_CHECKED_Minus].size() + PD_Minus_Output[TOTAL_SNP_ERROR_CHECKED_Minus].size()) {
				const short CurrentLengthOutput = CurrentLength + 1;
				CheckBoth(TheInput, CurrentReadSeq, PD_Plus_Output, PD_Minus_Output,
							 BP_Start, BP_End,
							 CurrentLengthOutput, UP);
			}
			else return;
		}
		 */
   }
   else return;
}


void CleanUniquePoints (vector <UniquePoint> & Input_UP) {
	vector <UniquePoint> TempUP; //vector <UniquePoint> UP_Close; UP_Far
	UniquePoint LastUP = Input_UP[Input_UP.size() - 1];
	//TempUP.push_back(LastUP);
	char LastDirection = LastUP.Direction;
	char LastStrand = LastUP.Strand;
	unsigned int Terminal;

	if (LastDirection == FORWARD) {
		Terminal = LastUP.AbsLoc - LastUP.LengthStr;
		for (unsigned i = 0; i < Input_UP.size(); i++) {
			if (Input_UP[i].Direction == LastDirection && Input_UP[i].Strand == LastStrand) {
				if (Terminal == Input_UP[i].AbsLoc - Input_UP[i].LengthStr)
					TempUP.push_back(Input_UP[i]);
			}
		}
	}
	else if (LastDirection == BACKWARD) {
		Terminal = LastUP.AbsLoc + LastUP.LengthStr;
		for (unsigned i = 0; i < Input_UP.size(); i++) {
			if (Input_UP[i].Direction == LastDirection && Input_UP[i].Strand == LastStrand) {
				if (Terminal == Input_UP[i].AbsLoc + Input_UP[i].LengthStr)
					TempUP.push_back(Input_UP[i]);
			}
		}
	}
   Input_UP.clear();
	Input_UP = TempUP;
}

void SortOutputLI(const string & CurrentChr, vector <SPLIT_READ> & Reads, ofstream & LargeInsertionOutf) {
	unsigned UP_Close_index;
	unsigned temp_AbsLoc;
	// find LI combinations
	uint8_t * plus_LI_Pos = new uint8_t[CurrentChr.size() + 1];
	uint8_t * minus_LI_Pos = new uint8_t[CurrentChr.size() + 1];
	for (unsigned i = 0; i < CurrentChr.size() + 1; i++) {
		plus_LI_Pos[i] = 0;
		minus_LI_Pos[i] = 0;
	}
	for (unsigned Index = 0; Index < Reads.size(); Index++) {
			UP_Close_index = Reads[Index].UP_Close.size() - 1;
			temp_AbsLoc = Reads[Index].UP_Close[UP_Close_index].AbsLoc;
			if (Reads[Index].MatchedD == Plus)
				plus_LI_Pos[temp_AbsLoc]++;
			else if (Reads[Index].MatchedD == Minus)
				minus_LI_Pos[temp_AbsLoc]++;
	}
	vector <LI_Pos> LI_Positions;
	LI_Pos temp_LI_pos;

	for (int Index_Plus = SpacerBeforeAfter; Index_Plus < CurrentChr.size() - SpacerBeforeAfter; Index_Plus++) {
		if (plus_LI_Pos[Index_Plus] >= NumRead2ReportCutOff) {
			for (int Index_Minus = Index_Plus - 30; Index_Minus <= Index_Plus + 1; Index_Minus++) {
				if (minus_LI_Pos[Index_Minus] >= NumRead2ReportCutOff) {
					temp_LI_pos.Plus_Pos = Index_Plus;
					temp_LI_pos.Minus_Pos = Index_Minus;
					temp_LI_pos.WhetherReport = false;
					LI_Positions.push_back(temp_LI_pos);
				}
         }
		}
	}
   //cout << "LI: " << LI_Positions.size() << endl;
	int Count_LI = 0;
	// find LI supporting reads
	for (unsigned Index = 0; Index < Reads.size(); Index++) {
		UP_Close_index = Reads[Index].UP_Close.size() - 1;
		temp_AbsLoc = Reads[Index].UP_Close[UP_Close_index].AbsLoc;
		for (unsigned LI_index = 0; LI_index < LI_Positions.size(); LI_index++) {
			if (Reads[Index].MatchedD == Plus) {
				if (temp_AbsLoc == LI_Positions[LI_index].Plus_Pos) {
					Reads[Index].Used = true;
					LI_Positions[LI_index].Plus_Reads.push_back(Index);
				}
			}
			else if (Reads[Index].MatchedD == Minus) {
				if (temp_AbsLoc == LI_Positions[LI_index].Minus_Pos) {
					Reads[Index].Used = true;
					LI_Positions[LI_index].Minus_Reads.push_back(Index);
				}
			}
		}
	}

	vector <SPLIT_READ> temp_Plus_Reads, temp_Minus_Reads;

	bool temp_BalancedPlus_Plus, temp_BalancedPlus_Minus, temp_BalancedMinus_Plus, temp_BalancedMinus_Minus;
	short temp_LengthStr;
	for (unsigned LI_index = 0; LI_index < LI_Positions.size(); LI_index++) {
		temp_BalancedPlus_Plus = false;
		temp_BalancedPlus_Minus = false;
		temp_BalancedMinus_Plus = false;
		temp_BalancedMinus_Minus = false;
		temp_Plus_Reads.clear();
		temp_Minus_Reads.clear();

		for (int i = 0; i < LI_Positions[LI_index].Minus_Reads.size(); i++) {
			temp_Minus_Reads.push_back(Reads[LI_Positions[LI_index].Minus_Reads[i]]);
			UP_Close_index = Reads[LI_Positions[LI_index].Minus_Reads[i]].UP_Close.size() - 1;
			temp_LengthStr = Reads[LI_Positions[LI_index].Minus_Reads[i]].UP_Close[UP_Close_index].LengthStr;
			if ((float)temp_LengthStr > Reads[LI_Positions[LI_index].Minus_Reads[i]].ReadLength * 0.5)
				temp_BalancedMinus_Plus = true;
			else if ((float)temp_LengthStr < Reads[LI_Positions[LI_index].Minus_Reads[i]].ReadLength * 0.5)
				temp_BalancedMinus_Minus = true;
		}
		for (int i = 0; i < LI_Positions[LI_index].Plus_Reads.size(); i++) {
			temp_Plus_Reads.push_back(Reads[LI_Positions[LI_index].Plus_Reads[i]]);
			UP_Close_index = Reads[LI_Positions[LI_index].Plus_Reads[i]].UP_Close.size() - 1;
			temp_LengthStr = Reads[LI_Positions[LI_index].Plus_Reads[i]].UP_Close[UP_Close_index].LengthStr;
			if ((float)temp_LengthStr > Reads[LI_Positions[LI_index].Plus_Reads[i]].ReadLength * 0.5)
				temp_BalancedPlus_Plus = true;
			else if ((float)temp_LengthStr < Reads[LI_Positions[LI_index].Plus_Reads[i]].ReadLength * 0.5)
				temp_BalancedPlus_Minus = true;
		}
		if (temp_BalancedPlus_Plus && temp_BalancedPlus_Minus && temp_BalancedMinus_Plus && temp_BalancedMinus_Minus) LI_Positions[LI_index].WhetherReport = true;
	}
	// output
	for (unsigned LI_index = 0; LI_index < LI_Positions.size(); LI_index++) {
		if (LI_Positions[LI_index].WhetherReport) {
			Count_LI++;
			LargeInsertionOutf << "########################################################" << endl;
			LargeInsertionOutf << "ChrID " << temp_Plus_Reads[0].FragName << "\t"
			<< LI_Positions[LI_index].Plus_Pos - SpacerBeforeAfter + 1 << "\t"
			<< temp_Plus_Reads.size() << "\t"
			<< LI_Positions[LI_index].Minus_Pos - SpacerBeforeAfter  + 1 << "\t"
			<< temp_Minus_Reads.size() << endl;

			LargeInsertionOutf <<  ( CurrentChr.substr(LI_Positions[LI_index].Plus_Pos - ReportLength + 1, ReportLength) )
			<< Cap2Low( CurrentChr.substr(LI_Positions[LI_index].Plus_Pos + 1, ReportLength ) ) << endl;
			for (int i = 0; i < temp_Plus_Reads.size(); i++) {
				UP_Close_index = temp_Plus_Reads[i].UP_Close.size() - 1;
				temp_LengthStr = temp_Plus_Reads[i].UP_Close[UP_Close_index].LengthStr;
				for (int j = 0; j < ReportLength - temp_LengthStr; j++) {
					LargeInsertionOutf << " ";
				}
				LargeInsertionOutf << ReverseComplement(temp_Plus_Reads[i].UnmatchedSeq) << endl;
			}

			LargeInsertionOutf << "--------------------------------------------------------" << endl;
			//LargeInsertionOutf << "-\t" << minus_LI_Pos[Index_Minus].NumReads << endl;
			LargeInsertionOutf << Cap2Low( CurrentChr.substr(LI_Positions[LI_index].Minus_Pos - ReportLength, ReportLength) )
			<< ( CurrentChr.substr(LI_Positions[LI_index].Minus_Pos, ReportLength ) ) << endl;
			for (int i = 0; i < temp_Minus_Reads.size(); i++) {
				UP_Close_index = temp_Minus_Reads[i].UP_Close.size() - 1;
				temp_LengthStr = temp_Minus_Reads[i].UP_Close[UP_Close_index].LengthStr;
				for (int j = 0; j < ReportLength + temp_LengthStr - temp_Minus_Reads[i].ReadLength; j++) {
					LargeInsertionOutf << " ";
				}
				LargeInsertionOutf << (temp_Minus_Reads[i].UnmatchedSeq) << endl;
			}
		}
	}

	cout << "Breakpoints for large insertions (LI): " << Count_LI << "\n\n";
}


void SortOutputRest(const string & CurrentChr, vector <SPLIT_READ> & Reads, ofstream & Outf_Rest) {
	unsigned UP_Close_index;
	unsigned temp_AbsLoc;
	// find LI combinations
	uint8_t * plus_LI_Pos = new uint8_t[CurrentChr.size() + 1];
	uint8_t * minus_LI_Pos = new uint8_t[CurrentChr.size() + 1];
	for (unsigned i = 0; i < CurrentChr.size() + 1; i++) {
		plus_LI_Pos[i] = 0;
		minus_LI_Pos[i] = 0;
	}
	for (unsigned Index = 0; Index < Reads.size(); Index++) {
		if (Reads[Index].Used) continue;
		UP_Close_index = Reads[Index].UP_Close.size() - 1;
		temp_AbsLoc = Reads[Index].UP_Close[UP_Close_index].AbsLoc;
		if (Reads[Index].MatchedD == Plus)
			plus_LI_Pos[temp_AbsLoc]++;
		else// if (Reads[Index].MatchedD == Minus)
			minus_LI_Pos[temp_AbsLoc]++;
	}
	vector <Rest_Pos> Rest_Positions;
	Rest_Pos temp_Rest_pos;

	for (int Index = SpacerBeforeAfter; Index < CurrentChr.size() - SpacerBeforeAfter; Index++) {
		if (plus_LI_Pos[Index] >= NumRead2ReportCutOff) {
			temp_Rest_pos.Strand = Plus;
			temp_Rest_pos.Pos = Index;
			Rest_Positions.push_back(temp_Rest_pos);
		}
		if (minus_LI_Pos[Index] >= NumRead2ReportCutOff) {
			temp_Rest_pos.Strand = Minus;
			temp_Rest_pos.Pos = Index;
			Rest_Positions.push_back(temp_Rest_pos);
		}
	}

	// find supporting reads
	for (unsigned Index = 0; Index < Reads.size(); Index++) {
		if (Reads[Index].Used) continue;
		UP_Close_index = Reads[Index].UP_Close.size() - 1;
		temp_AbsLoc = Reads[Index].UP_Close[UP_Close_index].AbsLoc;
		for (unsigned Pos_index = 0; Pos_index < Rest_Positions.size(); Pos_index++) {
			if (Reads[Index].MatchedD == Rest_Positions[Pos_index].Strand) {
				if (temp_AbsLoc == Rest_Positions[Pos_index].Pos) {
					Reads[Index].Used = true;
					Rest_Positions[Pos_index].Pos_Reads.push_back(Index);	 // copy index to save memory
				}
			}
		}
	}

	//cout << "Other unassigned breakpoints (BP): " << Rest_Positions.size() << "\n\n";
	int Count_BP = 0;
	bool temp_BalancedPlus, temp_BalancedMinus;
	short temp_LengthStr;
	vector <SPLIT_READ> temp_Pos_Reads;
	for (unsigned LI_index = 0; LI_index < Rest_Positions.size(); LI_index++) {
		temp_Pos_Reads.clear();
		temp_BalancedPlus = false;
		temp_BalancedMinus = false;
		for (int i = 0; i < Rest_Positions[LI_index].Pos_Reads.size(); i++) {
			temp_Pos_Reads.push_back(Reads[Rest_Positions[LI_index].Pos_Reads[i]]);
			UP_Close_index = Reads[Rest_Positions[LI_index].Pos_Reads[i]].UP_Close.size() - 1;
			temp_LengthStr = Reads[Rest_Positions[LI_index].Pos_Reads[i]].UP_Close[UP_Close_index].LengthStr;
			if ((float)temp_LengthStr > Reads[Rest_Positions[LI_index].Pos_Reads[i]].ReadLength * 0.5)
				temp_BalancedPlus = true;
			else if ((float)temp_LengthStr < Reads[Rest_Positions[LI_index].Pos_Reads[i]].ReadLength * 0.5)
				temp_BalancedMinus = true;
		}
		if (temp_BalancedPlus && temp_BalancedMinus) {
			Count_BP++;
			if (Rest_Positions[LI_index].Strand == Plus) {
				Outf_Rest << "########################################################" << endl;
				Outf_Rest << "ChrID " << temp_Pos_Reads[0].FragName << "\t"
				<< Rest_Positions[LI_index].Pos - SpacerBeforeAfter + 1 << "\t"
				<< Rest_Positions[LI_index].Pos_Reads.size() << endl;

				Outf_Rest <<  ( CurrentChr.substr(Rest_Positions[LI_index].Pos - ReportLength + 1, ReportLength) )
				<< Cap2Low( CurrentChr.substr(Rest_Positions[LI_index].Pos + 1, ReportLength ) ) << endl;
				for (int i = 0; i < temp_Pos_Reads.size(); i++) {
					UP_Close_index = temp_Pos_Reads[i].UP_Close.size() - 1;
					temp_LengthStr = temp_Pos_Reads[i].UP_Close[UP_Close_index].LengthStr;
					for (int j = 0; j < ReportLength - temp_LengthStr; j++) {
						Outf_Rest << " ";
					}
					Outf_Rest << ReverseComplement(temp_Pos_Reads[i].UnmatchedSeq)
					          << "\t" << temp_Pos_Reads[i].MatchedD
					          << "\t" << temp_Pos_Reads[i].MatchedRelPos
					          << "\t" << temp_Pos_Reads[i].MS
					          << "\t" << temp_Pos_Reads[i].Tag
					          << "\t" <<  temp_Pos_Reads[i].Name << endl;
				}

			}
			else {
				Outf_Rest << "########################################################" << endl;
				Outf_Rest << "ChrID " << temp_Pos_Reads[0].FragName << "\t"
				<< Rest_Positions[LI_index].Pos - SpacerBeforeAfter + 1 << "\t"
				<< Rest_Positions[LI_index].Pos_Reads.size() << endl;
				Outf_Rest << Cap2Low( CurrentChr.substr(Rest_Positions[LI_index].Pos - ReportLength, ReportLength) )
				<< ( CurrentChr.substr(Rest_Positions[LI_index].Pos, ReportLength ) ) << endl;
				for (int i = 0; i < temp_Pos_Reads.size(); i++) {
					UP_Close_index = temp_Pos_Reads[i].UP_Close.size() - 1;
					temp_LengthStr = temp_Pos_Reads[i].UP_Close[UP_Close_index].LengthStr;
					for (int j = 0; j < ReportLength + temp_LengthStr - temp_Pos_Reads[i].ReadLength; j++) {
						Outf_Rest << " ";
					}
					Outf_Rest << (temp_Pos_Reads[i].UnmatchedSeq)
					          << "\t" << temp_Pos_Reads[i].MatchedD
					          << "\t" << temp_Pos_Reads[i].MatchedRelPos
					          << "\t" << temp_Pos_Reads[i].MS
					          << "\t" << temp_Pos_Reads[i].Tag
					          << "\t" <<  temp_Pos_Reads[i].Name << endl;				}
			}
		}
	}

	cout << "Other unassigned breakpoints (BP): " << Count_BP << "\n\n";
}

short CompareTwoString(const string & Str_A, const string & Str_B) {
	short Str_Len = Str_A.size();
	short CompareResult;
	for (short i = 0; i < Str_Len; i++) {
		CompareResult = (short)Str_A[i] - (short)Str_B[i];
		if (CompareResult < 0) return -1;
		else if (CompareResult > 0) return 1;
	}
	return 0;
}

string GetConsensusInsertedStr(const vector <SPLIT_READ> & Reads, const int & StartIndex, const int & EndIndex) {
    // InsertedStr
    map<string,int> NT_str_2_count;
    map<string,int>::iterator it;
    for (int i = StartIndex; i <= EndIndex; i++) {
        it = NT_str_2_count.find(Reads[i].InsertedStr);
        if (it == NT_str_2_count.end()) {
            NT_str_2_count[Reads[i].InsertedStr] = 1;
        }
        else it->second++;
    }
    int Max = 0;
    string OutputStr = "";
    for (std::map<std::string,int>::iterator it=NT_str_2_count.begin(); it!=NT_str_2_count.end(); it++ ) {
        if (it->second > Max) {
            Max = it->second;
            OutputStr = it->first;
        }
    }
    return OutputStr;
}
