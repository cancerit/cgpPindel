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

short MAX_SNP_ERROR = 2;
short ADDITIONAL_MISMATCH = 2;

short TOTAL_SNP_ERROR_CHECKED = MAX_SNP_ERROR + ADDITIONAL_MISMATCH + 1;
short TOTAL_SNP_ERROR_CHECKED_Minus = MAX_SNP_ERROR + ADDITIONAL_MISMATCH;

short MAX_ALLOWED_MISMATCHES = TOTAL_SNP_ERROR_CHECKED_Minus + 5;

// #########################################################
short Min_Perfect_Match_Around_BP = 3;                   //#
const short MIN_IndelSize_NT = 3;                        //#
const short MIN_IndelSize_Inversion = 3;                 //#
float Seq_Error_Rate = 0.05;                             //#
float Seq_Error_Rate_1 = 0.05;                           //#
float Seq_Error_Rate_2 = 0.02;                           //#
float Seq_Error_Rate_3 = 0.00;                           //#
unsigned int BalanceCutoff = 50;                         //#
short RangeMaxSensivity = 9;       // 3                  //#
short RangeMediumSensivity = 9;    // 5                  //#
//short RangeLowSensivity = 7;                           //#
                                                         //#
const unsigned int NumRead2ReportCutOff = 3;             //#
short MaxRangeIndex = 9;// 5 or 6 or 7 or maximum 8//#
const float MaximumAllowedMismatchRate = 0.1;            //#
const short Min_Num_Matched_Bases = 30;                  //#
// #########################################################
//const float Double_Seq_Error_Rate_Per_Side = Seq_Error_Rate_Per_Side * 2;
unsigned int Distance = 300;
 short MinClose = 8;//short(log((double)Distance)/log(4.0) + 0.8) + 3 + MAX_SNP_ERROR;//atoi(argv[1]);
 short MinFar_I = MinClose + 1;//atoi(argv[2]);
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
short ReadInRead(ifstream & inf_Seq, const string & ChrSeq, const string & ChrName, ofstream & Outf);
void GetCloseEnd(const string & CurrentChr, SPLIT_READ & Temp_One_Read);
vector <string> ReverseComplement(const vector <string> & input);
string Reverse(const string & InputPattern);
string ReverseComplement(const string & InputPattern);
string Cap2Low(const string & input);
void CheckLeft_Close(const string & TheInput,
               const string & CurrentReadSeq,
               const vector <unsigned int> Left_PD[],
               const short & BP_Left_Start,
               const short & BP_Left_End,
               const short & CurrentLength,
               vector <UniquePoint> & LeftUP);

void CheckRight_Close(const string & TheInput,
                const string & CurrentReadSeq,
                const vector <unsigned int> Right_PD[],
                const short & BP_Right_Start,
                const short & BP_Right_End,
                const short & CurrentPos,
                vector <UniquePoint> & RightUP);

bool CheckMismatches(const string & TheInput,
                     const string & CurrentReadSeq,
                     //const unsigned int & Start,
							const UniquePoint & UP);

void CleanUniquePoints (vector <UniquePoint> & Input_UP);
int main (int argc, char *argv[]) {

	if (NumRead2ReportCutOff == 1) BalanceCutoff = 3000000000;

	// ################### module 1: define input and output #####################
  if (argc < 4) {
    cout << "\nFiltering Pindel reads, developed by Kai Ye, k.ye@lumc.nl\n\n"
         << "at least 4 parameters are required here:\n"
         << "1. Input: the reference genome sequences in fasta format;\n"
	      << "2. Which chr/fragment\n"
	      << "3. Output file name\n"
  	      << "4. Input(s): the unmapped reads in a modified fastq format;\n"
         << endl;
    return 1;
  }

  // #################################################################
	const string WhichChr = argv[2];

	ifstream  inf_Seq(argv[1]);   // input file name

	string OutputFile = argv[3];
	ofstream Output(OutputFile.c_str());
	if (!Output) {
		cout << "Sorry, cannot write to the file: " << OutputFile << endl;
		return 1;
	}


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

   string CurrentChr;

       ReadInOneChr(inf_Seq, CurrentChr, WhichChr);
	    if (CurrentChr.empty()) {
			 cout << "Cannot find the requested chr." << endl;
			 return 1;
		 }
	    CONS_Chr_Size = CurrentChr.size() - 2 * SpacerBeforeAfter;

	for (int FileIndex = 4; FileIndex < argc; FileIndex++) {
		cout << "processing file: " << argv[FileIndex] << endl;
		ifstream  inf_ReadsSeq(argv[FileIndex]);   // input file name
		ReadInRead(inf_ReadsSeq, CurrentChr, WhichChr, Output);
	}

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

short ReadInRead(ifstream & inf_ReadSeq, const string & CurrentChr, const string & FragName, ofstream & Outf) {
	cout << "Scanning and processing reads anchored in " << FragName << endl;
	//short ADDITIONAL_MISMATCH = 1;
	ADDITIONAL_MISMATCH = 1;
	Seq_Error_Rate = 0.05;
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
			//Reads.clear();
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
			MAX_SNP_ERROR = (short)(Temp_One_Read.UnmatchedSeq.size() * Seq_Error_Rate);

			TOTAL_SNP_ERROR_CHECKED = MAX_SNP_ERROR + ADDITIONAL_MISMATCH + 1;
			TOTAL_SNP_ERROR_CHECKED_Minus = MAX_SNP_ERROR + ADDITIONAL_MISMATCH;
			MinClose = short(log((double)(Temp_One_Read.InsertSize * 3))/log(4.0) + 0.8) + 3;// + MAX_SNP_ERROR;//atoi(argv[1]);
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
					//Reads.push_back(Temp_One_Read);
				Outf << Temp_One_Read.Name << "\n"
				<< Temp_One_Read.UnmatchedSeq << "\n"
				<< Temp_One_Read.MatchedD << "\t"
				<< Temp_One_Read.FragName << "\t"
				<< Temp_One_Read.MatchedRelPos << "\t"
				<< Temp_One_Read.MS << "\t"
				<< Temp_One_Read.InsertSize << "\t"
				<< Temp_One_Read.Tag << "\n";
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
		}
	}

	cout << "NumReadScanned:\t" << NumReadScanned << endl;
	cout << "NumReadInChr:\t" << NumReadInChr << endl;
	cout << "NumReadStored:\t" << GetPlus + GetMinus << endl;
	cout << "NumReadStored / NumReadInChr = " << (GetPlus + GetMinus) * 100.0 / NumReadInChr << " %\n"
	<< "InChrPlus \t" << InChrPlus << "\tGetPlus \t" << GetPlus << "\t" << GetPlus * 100.0 / InChrPlus
	<< " %\n" << "InChrMinus\t" << InChrMinus  << "\tGetMinus\t" << GetMinus
	<< "\t" << GetMinus * 100.0 / InChrMinus << " %\n" << endl;
	inf_ReadSeq.close();
	//if (Reads.size() == 0) return 0;
	//cout << LeftReads.size() << endl;
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

void CheckLeft_Close(const string & TheInput,
							const string & CurrentReadSeq,
							const vector <unsigned int> Left_PD[],
							const short & BP_Left_Start,
							const short & BP_Left_End,
							const short & CurrentLength,
							vector <UniquePoint> & LeftUP) {
	int Sum;
   if (CurrentLength >= BP_Left_Start && CurrentLength <= BP_Left_End) {
      // put it to LeftUP if unique
		for (short i = 0; i <= MAX_SNP_ERROR; i++) {
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
      vector <unsigned int> Left_PD_Output[TOTAL_SNP_ERROR_CHECKED];
		for (int CheckedIndex = 0; CheckedIndex < TOTAL_SNP_ERROR_CHECKED; CheckedIndex++) {
			Left_PD_Output[CheckedIndex].reserve(Left_PD[CheckedIndex].size());
		}
      const char CurrentChar = CurrentReadSeq[CurrentLength];
      //const int SizeOfCurrent = Left_PD.size();
      unsigned int pos;
		for (int i = 0; i < TOTAL_SNP_ERROR_CHECKED_Minus; i++) {
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
		Sum = 0;
		for (int i = 0; i <= MAX_SNP_ERROR; i++) {
			Sum += Left_PD_Output[i].size();
		}
      if (Sum) {
         const short CurrentLengthOutput = CurrentLength + 1;
         CheckLeft_Close(TheInput, CurrentReadSeq, Left_PD_Output,
								 BP_Left_Start, BP_Left_End,
								 CurrentLengthOutput, LeftUP);
      }
      else return;
   }
   else return;
}

void CheckRight_Close(const string & TheInput,
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
	   for (short i = 0; i <= MAX_SNP_ERROR; i++) {
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
      vector <unsigned int> Right_PD_Output[TOTAL_SNP_ERROR_CHECKED];
		for (int CheckedIndex = 0; CheckedIndex < TOTAL_SNP_ERROR_CHECKED; CheckedIndex++) {
			Right_PD_Output[CheckedIndex].reserve(Right_PD[CheckedIndex].size());
		}
      const char CurrentChar = CurrentReadSeq[ReadLengthMinus - CurrentLength];
      unsigned int pos;

		for (int i = 0; i < TOTAL_SNP_ERROR_CHECKED_Minus; i++) {
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

		Sum = 0;
		for (int i = 0; i <= MAX_SNP_ERROR; i++) {
			Sum += Right_PD_Output[i].size();
		}
      if (Sum) {
         short CurrentLength_output = CurrentLength + 1;
         CheckRight_Close(TheInput, CurrentReadSeq, Right_PD_Output,
								  BP_Right_Start, BP_Right_End,
								  CurrentLength_output, RightUP);
      }
      else return;
   }
   else return;
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

	MAX_ALLOWED_MISMATCHES = (short)(CurrentReadSeq.size() * MaximumAllowedMismatchRate + 1); //

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

void GetCloseEnd(const string & CurrentChr, SPLIT_READ & Temp_One_Read) {

	Temp_One_Read.ReadLength = Temp_One_Read.UnmatchedSeq.size();
	Temp_One_Read.ReadLengthMinus = Temp_One_Read.ReadLength - 1;
	char LeftChar, RightChar;
	string CurrentReadSeq;
	vector <unsigned int> PD[TOTAL_SNP_ERROR_CHECKED];
	for (int CheckIndex = 0; CheckIndex < TOTAL_SNP_ERROR_CHECKED; CheckIndex++) {
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
	MinClose = short(log((double)(Temp_One_Read.InsertSize * 3))/log(4.0) + 0.8) + 3;
	BP_Start = MinClose;
	BP_End = Temp_One_Read.ReadLengthMinus;
	//Temp_One_Read.OK = true;
	if (TOTAL_SNP_ERROR_CHECKED_Minus) {
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
			CheckLeft_Close(CurrentChr, CurrentReadSeq, PD, BP_Start, BP_End, FirstBase, UP);
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
			CheckRight_Close(CurrentChr, CurrentReadSeq, PD, BP_Start, BP_End, FirstBase, UP);
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
			CheckLeft_Close(CurrentChr, CurrentReadSeq, PD, BP_Start, BP_End, FirstBase, UP);
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
			CheckRight_Close(CurrentChr, CurrentReadSeq, PD, BP_Start, BP_End, FirstBase, UP);
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

