
/*  *** 2D CHAIN *** 
	Create Rigid Planar Chains on Honeycomb Lattice
	Check for Overlaps, Closed Loop Chains
	Systematic Approach, Checks All Chains of Given Length (less than 64 segments)
	GPU/CUDA Version, January 2014

	Rewritten to reflect number of segments rather than atoms, added minor improvements - March 2014

	By Christian Bracher */

#include "cuda_runtime.h"
#include "cuda_profiler_api.h"
#include "device_launch_parameters.h"
#include <stdio.h>
#include <iostream>

/* Support for Timing */

#include <time.h> 

/* Check for 64-bit support (necessary) */

#if defined (_INTEGRAL_MAX_BITS) && \
  _INTEGRAL_MAX_BITS >= 64
typedef signed __int64 int64;
typedef unsigned __int64 uint64;
#else
#error __int64 type not supported
#endif

using namespace std;

/* *** GPU CONTROL PARAMETERS *** */

/* Threads per Block (default: 512, maximum: 1024) */
	int ThreadBits = 9;
	int ThreadsPerBlock = (int)(1 << ThreadBits);
	const int MaxThreadsPerBlock = 1024;

/* Number of Blocks (default: 4096)  */
	int BlockBits = 12;
	long NumberOfBlocks = ((long)1 << BlockBits);
	const long MaxNumberOfBlocks = 65536;

	
/* *** LATTICE VECTOR CLASS AND FUNCTIONS *** */

class LatticeVector
{
/* Grid Coordinates on 2D Triangular Grid */
public:
	int n1;
	int n2;

	__device__ LatticeVector ()
	{
	/* Default Constructor */
		n1 = 0;
		n2 = 0;
	}

	__device__ LatticeVector (int x, int y)
	{
	/* Create a Lattice Point */
		n1 = x;
		n2 = y;
	}

	__device__ void LeftTurn (int orientation)
	{
	/* Calculates Grid Position of New Terminal Atom for a Left Turn, Assuming Bond Orientation d */
		
		int d = orientation % 6;
		if (d < 0) d += 6;

		switch (d)
		{
		case 0:
			++n2;
			break;
		case 1:
			--n1;
			++n2;
			break;
		case 2:
			--n1;
			break;
		case 3:
			--n2;
			break;
		case 4:
			++n1;
			--n2;
			break;
		case 5:
			++n1;
			break;
		}
	}

	__device__ void RightTurn (int orientation)
	{
	/* Calculates Grid Position of New Terminal Atom for a Right Turn, Assuming Bond Orientation d */
		
		int d = orientation % 6;
		if (d < 0) d += 6;

		switch (d)
		{
		case 0:
			++n1;
			--n2;
			break;
		case 1:
			++n1;
			break;
		case 2:
			++n2;
			break;
		case 3:
			--n1;
			++n2;
			break;
		case 4:
			--n1;
			break;
		case 5:
			--n2;
			break;
		}
	}

	__device__ bool operator == (LatticeVector r1)
	{
	/* Compares Two Lattice Vectors */

		bool IsEqual = false;

		if ((r1.n1 == n1) && (r1.n2 == n2))  (IsEqual = true);

		return IsEqual;
	}
};

__device__ void BuildChain (uint64 Code, int length, LatticeVector *ChainArray)
{
/* Translate the Binary Code Into the Actual Lattice Points Occupied by the Chain */

	ChainArray[0] = LatticeVector(0,0);
	ChainArray[1] = LatticeVector(1,0);
	ChainArray[2] = LatticeVector(1,1);

	int orientation = 1;

	for (int k = 2; k < length; ++k)
	{
		ChainArray[k+1] = ChainArray[k];

		if (((Code>>(k-2)) % 2) == 0)
		{
			ChainArray[k+1].LeftTurn(orientation);
			++orientation;
		}
		else
		{
			ChainArray[k+1].RightTurn(orientation);
			--orientation;
		}
	}
}

__device__ bool IsChainOverlapping (int length, LatticeVector *ChainArray)
{
/* Compare All Pairs of "Atoms" for Occupying the Same Grid Position
   (note that it is impossible to form loops with less than six atoms,
   overlapping vertices are always an even number of segments apart) */

	for (int k1 = 0; k1 < length - 5 ; ++k1)
	{
		for (int k2 = k1 + 6; k2 <= length; k2 += 2)
		{
			if (ChainArray[k2]==ChainArray[k1]) return true;
		}	
	}

	return false;
}

__device__ bool IsChainClosedLoop (int length, LatticeVector *ChainArray)
{
/* For a Closed Loop, No Atoms Overlap Except the First and Last */

	if (!(ChainArray[length] == ChainArray[0]))
	{
	/* Check Ends of Chain First - Are They Different? */
		return false;
	}
	else
	{
	/* Check for Overlaps in the Interior of the Chain */
		for (int k1 = 0; k1 < length - 5 ; ++k1)
		{
			for (int k2 = k1 + 6; k2 <= length; k2 += 2)
			{
				if (ChainArray[k2]==ChainArray[k1])
				{
					if ((k1 > 0) || (k2 < length))
						return false;
				}
			}	
		}
	}
	return true;
}

/* *** The ChainCounts Structure *** */
	
struct ChainCounts
	{
		unsigned long NonOverlapping;
		unsigned long ClosedChains;
	};

/* NEEDS TO BE ADAPTED... */ 

void PrintChainArray (int length, LatticeVector* ChainArray)
{
/* Print the Coordinates of the Atomic Grid Sites in the Chain */

	for(int i=0; i <= length; ++i)
	{
		cout << "(" << ChainArray[i].n1 << "," << ChainArray[i].n2 << ") ";
	}
	cout << "\n";
}

/* *** Parallel Code for Chain Analysis *** */

__global__ void CUDAChainAnalyze(ChainCounts *CUDAChainInfo, uint64 CodeOffset, int ChainLength)
{
/* Prepare Block Cache For Overlap Data */
	__shared__ ChainCounts DataCache[MaxThreadsPerBlock];

/* Allocate Memory for Building Chains (maximum length: 64) */
	LatticeVector MyChainArray[64];

/* Figure Out Correct Code For Two-Dimensional Chain */
	uint64 MyCode = CodeOffset + (uint64)(threadIdx.x + blockIdx.x * blockDim.x);	

/* Build the Chain ... */
	BuildChain(MyCode, ChainLength, MyChainArray);

/* ... and Test It For Overlaps (default value: overlaps): */
	DataCache[threadIdx.x].NonOverlapping = 0;
	
	if (IsChainOverlapping(ChainLength, MyChainArray) == false)
		DataCache[threadIdx.x].NonOverlapping = 1;
		
/* ... and Test It For Closed Chains (default value: not closed): */
	DataCache[threadIdx.x].ClosedChains = 0;
	
	if (IsChainClosedLoop(ChainLength, MyChainArray) == true)
		DataCache[threadIdx.x].ClosedChains = 1;
		
/* Wait for Tests Within a Block To Be Completed */
	__syncthreads();

/* Now, Add Results Within Block */
	int AddLimit = blockDim.x / 2;
	
	while (AddLimit > 0)
	{
		if (threadIdx.x < AddLimit)
		{
			DataCache[threadIdx.x].NonOverlapping += DataCache[threadIdx.x + AddLimit].NonOverlapping;
			DataCache[threadIdx.x].ClosedChains   += DataCache[threadIdx.x + AddLimit].ClosedChains;
		}

		__syncthreads();

		AddLimit /= 2;
	}

/* Store Result (now in position 0) in Global Memory */
	if (threadIdx.x == 0)
	{
		CUDAChainInfo[blockIdx.x].NonOverlapping = DataCache[0].NonOverlapping;
		CUDAChainInfo[blockIdx.x].ClosedChains = DataCache[0].ClosedChains;
	}
}

/* *** TIMING FUNCTONS *** */

double Duration (clock_t initial, clock_t final)
{
/* Determine the Duration of a Calculation in Seconds */

	double TimeDiff = (double)(final - initial);
	return TimeDiff / CLOCKS_PER_SEC;
}

/* *** MAIN PROGRAM STARTS HERE *** */

int main()
{
	int Length;
	clock_t StartTime, FinishTime;
	double StartToFinish;

/* Enter Chain Length Examined */
	cout << "Enter Chain Length: ";
	cin >> Length;
	cout << "\n\n";

/* Timing Support - Start of Calculation */
	StartTime = clock();

//	cudaProfilerStart();
		
/* Message */
	cout << "Calculating chains of length " << Length << " : ";

/* There are 2^(l-2) Different Chains */ 
	uint64 MaxCode = ((uint64)1 << (Length - 2));

/* Chop GPU Calculation Into Pieces If Length Exceeds ThreadBits + BlockBits + 2 */
	uint64 CUDAIterations;
	
	if (MaxCode > (NumberOfBlocks * ThreadsPerBlock))
	{
	/* Number of Global Iterations */
		CUDAIterations = ((MaxCode >> BlockBits) >> ThreadBits);
	}
	else
	{
	/* Adjust Number of Blocks */
		CUDAIterations = 1;
		NumberOfBlocks = (long)(MaxCode >> ThreadBits);

	/* Check: There Must Be At Least As Many Chains As Threads */
		if (NumberOfBlocks == 0)
		{
			cerr << "ERROR:  Insufficient Length of Chain \n\n";
			exit(1);
		}
	}
	
/* Initialize Variable to Save Numbers of Non-Overlapping Chains, Closed Chains */
	uint64 NonOverlapping = 0;
	uint64 ClosedChains = 0;

/* Reserve Memory to Transfer Information Between CPU and GPU */
	ChainCounts *CPUChainInfo = new ChainCounts[NumberOfBlocks];

	ChainCounts *CUDAChainInfo; 
	cudaError_t cudaStatus = cudaMalloc((void**)&CUDAChainInfo, NumberOfBlocks * sizeof(ChainCounts));

/* Loop Through All Possible Configurations */
	for (uint64 Iter = 0; Iter < CUDAIterations; ++Iter)
	{
	/* Examine (ThreadsPerBlock * BlockNumber) Chains in Parallel */

	/* Determine Offset for Parallel Calculation */
		uint64 Offset = Iter * NumberOfBlocks * ThreadsPerBlock;

	/* Perform Parallel Analysis */
		CUDAChainAnalyze <<<NumberOfBlocks,ThreadsPerBlock>>> (CUDAChainInfo, Offset, Length);

	/* Copy Results to Host Memory */
		cudaStatus = cudaMemcpy(CPUChainInfo, CUDAChainInfo, NumberOfBlocks * sizeof(ChainCounts), cudaMemcpyDeviceToHost);

	/* Extract Information */
		for (uint64 BlockID = 0; BlockID < NumberOfBlocks; ++BlockID)
		{
			NonOverlapping += (uint64)CPUChainInfo[BlockID].NonOverlapping;
			ClosedChains   += (uint64)CPUChainInfo[BlockID].ClosedChains; 
		}

	/* Indicate Progress */
		cout << ".";
	}

	cudaDeviceReset();

/* Timing Support - End of Calculation */
	FinishTime = clock();
	StartToFinish = Duration(StartTime, FinishTime);

/* Send a Brief Message */
		cout << " done!\n\n"
			 << "Found " << NonOverlapping << " non-overlapping chains.\n"
			 << "Found " << ClosedChains << " closed chains \n\n"
			 << "Time of Calculation: " << StartToFinish << " seconds.\n\n";

/* Wait for Key: */
	char Aux;
	cout << "Hit A Key, Then ENTER\n\n";
	cin >> Aux; 

/* Cleanup & Done! */
	delete[] CPUChainInfo;
	cudaFree(CUDAChainInfo);

	return 0;
}