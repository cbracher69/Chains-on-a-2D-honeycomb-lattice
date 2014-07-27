/*  *** 2D CHAIN ***
	Create Rigid Planar Chains on Honeycomb Lattice
	Check for Overlaps, Closed Loop Chains
	Systematic Approach, Checks All Chains of Given Length (less than 64 segments)

	Version 1:
	Two-dimensional grid approach, brute force check for overlaps in chain.

	Version 2.0:
	Implement Embedding Model in Three Dimensions

	Version 2.1:
	Include Smarter Overlap Check by Distance

	Version 2.2:
	Lower Bits/Digits in Code Now Represent the Free End of the Chain

	Version 2.3:
	Allow Building of Partial Chains - Reconstruct Free Ends

	Version 2.4:
	Smart Skip Algorithm for Overlapping Chains

	Version 2.5:
	Avoid Unnecessary Checks for Overlaps

	Version 3.1:
	Enable Self-Avoiding Polygon Codes

	Version 3.2:
	Sort Polygons, Eliminate Duplicates

	Version 3.3, 3.4:
	Examine Symmetry Classes of Polygons

	Version 3.5, 3.6:
	Improved Bound for Number of Closed Polygons

	CPU Version, January 2014 - May 2014
	By Christian Bracher */

#include "stdafx.h"
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

class LatticeVector
{
/* Grid Coordinates on 3D Cubic Grid, Used For Embedding Honeycomb Lattice */
public:
	int n1;
	int n2;
	int n3;

	LatticeVector ()
	{
	/* Default Constructor */
		n1 = 0;
		n2 = 0;
		n3 = 0;
	}

	LatticeVector (int x, int y, int z)
	{
	/* Create a Lattice Point */
		n1 = x;
		n2 = y;
		n3 = z;
	}

	inline void LeftTurn (int orientation)
	{
	/* Calculates Grid Position of New Terminal Atom for a Left Turn, Assuming Bond Orientation d */

		int d = orientation % 6;
		if (d < 0) d += 6;

		switch (d)
		{
		case 0:
			--n3;
			break;
		case 1:
			++n2;
			break;
		case 2:
			--n1;
			break;
		case 3:
			++n3;
			break;
		case 4:
			--n2;
			break;
		case 5:
			++n1;
			break;
		}
	}

	inline void RightTurn (int orientation)
	{
	/* Calculates Grid Position of New Terminal Atom for a Right Turn, Assuming Bond Orientation d */

		int d = orientation % 6;
		if (d < 0) d += 6;

		switch (d)
		{
		case 0:
			--n2;
			break;
		case 1:
			++n1;
			break;
		case 2:
			--n3;
			break;
		case 3:
			++n2;
			break;
		case 4:
			--n1;
			break;
		case 5:
			++n3;
			break;
		}
	}

	inline bool operator == (LatticeVector r1)
	{
	/* Compares Two Lattice Vectors */

		if ((r1.n1 == n1) && (r1.n2 == n2) && (r1.n3 == n3))
			return true;

		return false;
	}

	inline int distance (LatticeVector r1)
	{
	/* Finds Distance (L1 norm) to Another Point in the Lattice */

		return abs(r1.n1 - n1) + abs(r1.n2 - n2) + abs(r1.n3 - n3);
	}
};

class PolyMath
{
/* Functions Manipulating Codes Representing Closed-Loop Chains */
public:
	uint64 Code;
	int Length;

/* Constructor */
	PolyMath ()
	{
	/* Default Constructor */
		Code = 0;
		Length = 0;
	}

	PolyMath (uint64 chain, int length)
	{
	/* Create a Complete Polygon Code from an Open Chain Code */

	/* First, Copy Number of Segments */
		Length = length;

	/* Implicit Starting Turn is Left Turn (Zero Bit), Make Space For Final Turn */
		Code = (chain << 1);

	/* Figure Out Final Turn (Reproducing Polygon) */

	/* Trace the Orientation of Segments Along the Chain */
	/* Orientation of Second Segment */
		int orientation = 1;

	/* Prepare Comparison Code for Segment Analysis */
		uint64 PositionCode = ((uint64)1 << (length - 2));

		for (int k = 3; k <= length; ++k)
		{
		/* Select Next Bit for Comparison */
			PositionCode >>= 1;

		/* Select Segment, Calculate Next Lattice Position */
			if ((chain & PositionCode) == 0)
			{
				++orientation;
			}
			else
			{
				--orientation;
			}
		}

	/* Analyze Orientation of Last Segment of Chain */

	/* Normalize Orientation First */
		int d = (orientation % 6);
		if (d < 0) d += 6;

	/* Right Turn? */
		if (d == 1)
		{
			++Code;
		}
		else
		{
			if  (d != 5)  cerr << "ERROR: " << d << " ";
		}

	}

/* Rotate Closed Chains By One Segment */

	void operator ++ (void)
	{
	/* Rotate Chain to the Right */

	/* Isolate Last Segment in Chain */
		int FinalStep = (Code & 1);

	/* Shift Chain One Segment To Right */
		Code >>= 1;

	/* Make Last Segment First Segment */
		if (FinalStep == 1)
		{
			Code |= ((uint64)1 << (Length - 1));
		}
	}

	void operator -- (void)
	{
	/* Rotate Chain to the Left */

	/* Isolate First Segment in Chain */
		uint64 FirstStep = (Code >> (Length - 1));

	/* Shift Chain One Segment To Left */
		Code <<= 1;

	/* Make First Segment Last Segment */
		if (FirstStep == 1)
		{
		/* Eliminate Overflow Bit in Code */
			Code %= ((uint64)1 << Length);

		/* Replace First Segment */
			++Code;
		}
	}

/* Revert A Polygon Chain */
/* (Note: Reversion = run through chain backwards & exchange left and right turns) */

	void Revert (void)
	{
	/* Revert Polygon Chain */
		uint64 RevertedCode = 0;

	/* Trace Through Segments (Beginning Toward End) */
		for (int k = 0; k < Length; ++k)
		{

		/* Find Source Segment */
			uint64 SourceSegment = (Code & ((uint64)1 << k));

		/* Add As Target Segment in Reverse */
			if (SourceSegment == 0)
			{
				RevertedCode |= ((uint64)1 << (Length - k - 1));
			}
		}

		Code = RevertedCode;
	}

/* Reflect A Polygon Chain */
/* (Create a Mirror Image of the Polygon) */

	void Reflect (void)
	{
	/* Create Mirror Image of Polygon (by exchanging left and right turns) */

	/* Create Template */
		uint64 InvertCode = ((uint64)1 << Length);
		--InvertCode;

	/* Exchange Turns via XOR Operation */
		Code ^= InvertCode;
	}

/* Find Primitive Code of Closed Chain (rotated/reverted version with smallest code number) */

	void Reduce (void)
	{
	/* Storage for Smallest Code Found */
		uint64 MinCode = (((uint64)1 << Length) - 1);

	/* *** Step #1: Build Reverted Chain Code (reversed reading direction, reversed turns) */

	/* Create Variable for Reverted Code */
		uint64 RevertedCode = 0;

	/* Trace Through Segments (Beginning Toward End) */
		for (int k = 0; k < Length; ++k)
		{

		/* Find Source Segment */
			uint64 SourceSegment = (Code & ((uint64)1 << k));

		/* Add As Target Segment in Reverse */
			if (SourceSegment == 0)
			{
				RevertedCode |= ((uint64)1 << (Length - k - 1));
			}
		}

	/* *** Step 2: Loop Through Rotated Versions of Chain, Reverted Chain */

		for (int j = 0; j < 2; ++j)
		{
		/* Select First Original, Then Reverted Chain */
			if (j > 0)  Code = RevertedCode;

			for (int i = 0; i < Length; ++i)
			{
			/* Check Codes */
				if (Code < MinCode)  MinCode = Code;

			/* Rotate Chain One Segment:  Isolate Last Segment in Chain */
				int FinalStep = (Code & 1);

			/* Shift Chain One Segment To Right */
				Code >>= 1;

			/* Make Last Segment First Segment */
				if (FinalStep == 1)
				{
					Code |= ((uint64)1 << (Length - 1));
				}
			}
		}

	/* Replace Code by Primitive Code */
		Code = MinCode;
	}

/* Symmetry Functions */

/* Find Rotational Symmetry Class */

	int RotationalSymmetry (void)
	{
	/* Create Copy of Code */
		uint64 WorkCode = Code;

	/* Rotational Symmetry Parameter */
		int PtGrp = 0;

	/* Loop Thru Rotated Versions of Chain */
		for (int i = 0; i < Length; ++i)
		{
		/* Check Code - Got Original Back? */
			if (Code == WorkCode)  ++PtGrp;

		/* Rotate Chain One Segment:  Isolate Last Segment in Chain */
			int FinalStep = (WorkCode & 1);

		/* Shift Chain One Segment To Right */
			WorkCode >>= 1;

		/* Make Last Segment First Segment */
			if (FinalStep == 1)
			{
				WorkCode |= ((uint64)1 << (Length - 1));
			}
		}

	/* Return Result */
		return PtGrp;
	}

/* Examine Mirror Symmetry */
	bool MirrorSymmetry (void)
	{
	/* *** Step #1: Build Inverted Chain Code (reversed reading direction) */

	/* Create Variable for Inverted Code */
		uint64 InvertCode = 0;

	/* Trace Through Segments (Beginning Toward End) */
		for (int k = 0; k < Length; ++k)
		{

		/* Find Source Segment */
			uint64 SourceSegment = (Code & ((uint64)1 << k));

		/* Add As Target Segment, If Necessary */
			if (SourceSegment > 0)
			{
				InvertCode |= ((uint64)1 << (Length - k - 1));
			}
		}

	/* *** Step #2:  Rotate & Compare */

	/* Flag for Symmetry */
		bool IsMirrorSymmetric = false;

	/* Loop Thru Rotated Versions of Chain */
		for (int i = 0; i < Length; ++i)
		{
		/* Check Code - Is a Rotated Version of the Inverted Code the Original Code? */
			if ((InvertCode == Code))
			{
				IsMirrorSymmetric = true;
				break;
			}

		/* Rotate Chain One Segment:  Isolate Last Segment in Chain */
			int FinalStep = (InvertCode & 1);

		/* Shift Chain One Segment To Right */
			InvertCode >>= 1;

		/* Make Last Segment First Segment */
			if (FinalStep == 1)
			{
				InvertCode |= ((uint64)1 << (Length - 1));
			}
		}

	/* Return Result */
		return IsMirrorSymmetric;
	}
};


void BuildChain (uint64 Code, int length, LatticeVector *ChainArray)
{
/* Build a Complete Chain
   Translate the Binary Code Into the Actual Lattice Points Occupied by the Chain
   "0" Indicates Left Turn, "1" Indicates Right Turn
   Lowest Bits Indicate Free End of Chain ("going backwards")
   length Denotes Now Total Number of Segments */

/* Let Every Chain Start From the Origin to the Right, Then Make a Left Turn */
	ChainArray[0] = LatticeVector(0,0,0);
	ChainArray[1] = LatticeVector(1,0,0);
	ChainArray[2] = LatticeVector(1,0,-1);

/* Orientation of Second Segment */
	int orientation = 1;

/* Prepare Comparison Code for Segment Analysis */
	uint64 PositionCode = ((uint64)1 << (length - 2));

	for (int k = 3; k <= length; ++k)
	{
	/* Copy Current Lattice Position */
		ChainArray[k] = ChainArray[k-1];

	/* Select Next Bit for Comparison */
		PositionCode >>= 1;

	/* Select Segment, Calculate Next Lattice Position */
		if ((Code & PositionCode) == 0)
		{
			ChainArray[k].LeftTurn(orientation);
			++orientation;
		}
		else
		{
			ChainArray[k].RightTurn(orientation);
			--orientation;
		}
	}
}

inline int BranchingSegment(uint64 Code1, uint64 Code2, int length)
{
/* Find the Position of the First Segment That Deviates Between Two Chains */

/* Find Difference in Bit Structure Using XOR */
	uint64 DifferenceMap = (Code1 ^ Code2);

/* Determine Starting Position of "Tail" To Be Changed, Using Bit Shifts */
	int StartPos = length + 1;

/* Find Largest Altered Bit */
	do
	{
		DifferenceMap >>= 1;
		--StartPos;
	}
	while (DifferenceMap > 0);

	return StartPos;
}


inline void RebuildChain (uint64 Code, int StartPos, int length, LatticeVector *ChainArray)
{
/* Reconstruct The Free End of an Existing Chain Under a Change of Code
   length Denotes Total Number of Segments
   StartPos Denotes the First Segment to Be Rebuilt */

/* *** Task #1:  Find Orientation of Segment at Start of Rebuild Section */

/* Orientation of Second Segment */
	int orientation = 1;

/* Prepare Comparison Code for Segment Analysis */
	uint64 PositionCode = ((uint64)1 << (length - 2));

	for (int k = 3; k < StartPos; ++k)
	{
	/* Select Next Bit for Comparison */
		PositionCode >>= 1;

	/* Calculate Orientation of Next Segment */
		if ((Code & PositionCode) == 0)
			++orientation;
		else
			--orientation;
	}

/* *** Task #2:  Reconstruct End of Chain */

	for (int k = StartPos; k <= length; ++k)
	{
	/* Copy Current Lattice Position */
		ChainArray[k] = ChainArray[k-1];

	/* Select Next Bit for Comparison */
		PositionCode >>= 1;

	/* Select Segment, Calculate Next Lattice Position */
		if ((Code & PositionCode) == 0)
		{
			ChainArray[k].LeftTurn(orientation);
			++orientation;
		}
		else
		{
			ChainArray[k].RightTurn(orientation);
			--orientation;
		}
	}
}

inline int ChainOverlap (int segment, int length, LatticeVector *ChainArray)
{
/* Find the First Overlap of Two "Atoms" in Chain of Length (length),
   Assuming There Is No Overlap in the Initial Part of the Chain Up to Atom# (segment)

   Function Value Returned is the Position of the Overlapping Atom
   Function Returns Zero Result If Chain Has No Overlaps

   Method: Compare Pairs of "Atoms" in Ascending Order, Check for Occupying the Same Grid Position

   Notes:
   It is impossible to form loops with less than six atoms.
   The number of segments between two lattice points is at least their L1 distance. */

/* Loop Through "Target Atoms" */
	for (int k1 = segment; k1 <= length ; ++k1)
	{
	/* Check Ascending Chain for Overlap of Target Atom */
		int k2 = 0;

		while (k2 < k1 - 5)
		{
		/* Find Distance in L1 Metric */
			int separation = ChainArray[k1].distance(ChainArray[k2]);

			if (separation == 0)
			{
			/* Found Overlap; Report Position */
				return k1;
			}
			else
			{
				k2 += separation;
			}
		}
	}

/* No Overlap Found */
	return 0;
}

bool ClosedLoopCheck (int length, LatticeVector *ChainArray)
{
/* A Simplified Closed Loop Check that Assumes that The Overlapping Atom is the Final Atom in the Chain */

/* Check for Overlaps in Interior of Chain */
	for (int k1 = length - 6; k1 > 0; --k1)
	{
		if (ChainArray[length] == ChainArray[k1]) return false;
	}

/* None Found, So the Loop Must Be Closed */
	return true;
}

bool IsChainClosedLoop (int length, LatticeVector *ChainArray)
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

void PrintChainArray (int length, LatticeVector* ChainArray)
{
/* Print the Coordinates of the Atomic Grid Sites in the Chain */

	for(int i=0; i <= length; ++i){
		cout << "(" << ChainArray[i].n1 << "," << ChainArray[i].n2 << "," << ChainArray[i].n3 << ") ";
	}
	cout << "\n";
}

double Duration (clock_t initial, clock_t final)
{
/* Determine the Duration of a Calculation in Seconds */

	double TimeDiff = (double)(final - initial);
	return TimeDiff / CLOCKS_PER_SEC;
}

int _tmain(int argc, _TCHAR* argv[])
{
	int Length;
	clock_t StartTime, FinishTime;
	double StartToFinish;

/* Enter Maximum Chain Length Examined */
	cout << "Enter Chain Length: ";
	cin >> Length;
	cout << "\n\n";

/* Storage Array for Atom Coordinates */
	LatticeVector *ChainArray = new LatticeVector [Length + 1];

/* Storage Array for Self-Avoiding Polygon Chains */
/* (Estimate of number from exponential fit of known results) */
	uint64 MaxPolygonNumber = 32 + (uint64)(0.005 * exp(0.574 * Length));
	PolyMath *Polygon = new PolyMath [MaxPolygonNumber];

/* There are 2^(l-2) Different Chains */
	uint64 MaxCode = ((uint64)1 << (Length - 2));

	int OverlapAt, Segment;
	uint64 NonOverlaps = 0;
	uint64 ClosedChains = 0;

/* Auxiliary Variables for Progress Report */
	uint64 ChainChecks = 0;
	uint64 ProgressMark = (uint64)(1 << 24) - 1;

/* EXAMINE ALL CHAINS FOR OVERLAPS */

/* Timing Support - Start of Calculation */
	cout << "Calculating chains of length " << Length << " ... ";
	StartTime = clock();

/* Build Initial Chain */
	uint64 Code = 0;
	BuildChain(Code, Length, ChainArray);

/* Use Final Chain for Initial Comparison
   (guarantees complete check of initial chain for overlaps) */
	uint64 LastCode = MaxCode - 1;

/* Loop Through Chains, Search for Overlaps and Closed Self-Avoiding Chains */
	do
	{
	/* *** Task #1:  Progress Indicator, Check Counter */

	/* Count Examined Chains: */
		++ChainChecks;

	/* Progress Report after 2^24 (about 16 million) Evaluations */
		if ((ChainChecks & ProgressMark) == 0)
		{
			cerr << 100 * ((float)Code)/MaxCode << "% done.\n";
		}

	/* *** Task #2:  Find Common Head of Old and New Chains, Rebuild Chain, Check for Overlaps */

	/* Find Branching Segment */
		Segment = BranchingSegment(Code, LastCode, Length);

	/* Reconstruct Chain As Necessary */
		RebuildChain(Code, Segment, Length, ChainArray);

	/* Test Chain For Overlaps: */
		OverlapAt = ChainOverlap(Segment, Length, ChainArray);

	/* *** Task #3:  Intelligent Update of Chain Code, Jump Over Known "Bad" Chains */

	/* Remember Code Used */
		LastCode = Code;

	/* Analyze Result: No Overlap? */
		if (OverlapAt == 0)
		{
			++NonOverlaps;
			++Code;
		}
		else
		{
		/* Analyze Result: Closed Non-Overlapping Chain? */
			if (OverlapAt == Length)
			{
				if (ClosedLoopCheck(Length, ChainArray) == true)
				{
					Polygon[ClosedChains] = PolyMath(Code, Length);
					++ClosedChains;
				}
			}

		/* Perform "Smart Jump" to Next Code Without This Overlap */

		/* Remove Trailing End of Code */
			Code >>= (Length - OverlapAt);
		/* Step to Next Chain Segment Without This Specific Overlap */
			++Code;
		/* Fill Chain with Left Turns:  Lowest Value of Code Without Offending Overlap */
			Code <<= (Length - OverlapAt);
		}
	}
	while (Code < MaxCode);

/* Timing Support - End of Calculation */
	FinishTime = clock();
	StartToFinish = Duration(StartTime, FinishTime);

/* Send a Brief Message */
	cout << " done! \n\n";
	cout << "(Evaluations performed: " << ChainChecks << " out of " << MaxCode << " in " << StartToFinish << " seconds) \n\n";

/* SORT OUT SELF-AVOIDING POLYGONS */

/* Timing Support - Start of Calculation */
	cout << "Now Examining " << ClosedChains << " Self-Avoiding Polygons ... \n\n";
	StartTime = clock();

/* *** Step #1:  Find & Replace Reduced Polygon Codes */
	cout << "Reduce to Primitives ... ";

	for (uint64 i = 0; i < ClosedChains; ++i)
	{
		Polygon[i].Reduce();
	}

	cout << "done.\n";

/* *** Step #2:  Sort List of Primitives (using a fast merge-sort algorithm) */
	cout << "Sort List of Primitives ... ";

/* Create Array for Polygon Codes */
	uint64 *CodeArray = new uint64[ClosedChains];

/* Copy Polygon Codes Into Code Array */
	for (uint64 i = 0; i < ClosedChains; ++i)
	{
		CodeArray[i] = Polygon[i].Code;
	}

/* Get Rid of the Original Polygon Array */
	delete[] Polygon;

/* Create Temporary Array for Sort */
	uint64 *WorkArray = new uint64[ClosedChains];

/* Outer Sort Loop:  Double Sub-Arrays To Be Sorted by Factor of 2 in Each Run*/
	for (uint64 ArraySize = 1; ArraySize <= ClosedChains; ArraySize <<= 1)
	{
	/* Select Adjacent Sub-Arrays for Comparison Sort */
		for (uint64 ArrayCount = 0; ArrayCount < ClosedChains; ArrayCount += 2 * ArraySize)
		{
		/* Compare Primitives in Sub-Arrays */

		/* Bounds in Work Array */
			uint64 WorkCtr = ArrayCount;
			uint64 WorkCtrMax = ArrayCount + 2 * ArraySize;
			if (WorkCtrMax > ClosedChains)  WorkCtrMax = ClosedChains;

		/* Starting Pointers in Code Array */
			uint64 ArrayCtr1 = ArrayCount;
			uint64 ArrayCtr1Max = ArrayCount + ArraySize;
			if (ArrayCtr1Max > ClosedChains)  ArrayCtr1Max = ClosedChains;

			uint64 ArrayCtr2 = ArrayCtr1Max;

		/* Do the Sort */
			do
			{
				if (ArrayCtr1 < ArrayCtr1Max)
				{
				/* Sub-Array #1 Has Elements Available */
					uint64 Code1 = CodeArray[ArrayCtr1];

					if (ArrayCtr2 < WorkCtrMax)
					{
					/* Sub-Array #2 Has Elements Available */
						uint64 Code2 = CodeArray[ArrayCtr2];

						if (Code1 < Code2)
						{
						/* Comparison of Sub-Array Elements */
							WorkArray[WorkCtr] = Code1;

						/* Increment Counter in First Sub-Array */
							++ArrayCtr1;
						}
						else
						{
						/* Select Element From Second Sub-Array */
							WorkArray[WorkCtr] = Code2;

						/* Increment Counter in Second Sub-Array */
							++ArrayCtr2;
						}
					}
					else
					{
					/* Select Element From First Sub-Array By Default */
						WorkArray[WorkCtr] = Code1;

					/* Increment Counter in First Sub-Array */
						++ArrayCtr1;
					}
				}
				else
				{
				/* Select Element From Second Sub-Array By Default */
					WorkArray[WorkCtr] = CodeArray[ArrayCtr2];

				/* Increment Counter in Second Sub-Array */
					++ArrayCtr2;
				}

			/* Increase Work Counter */
				++WorkCtr;
			}
			while (WorkCtr < WorkCtrMax);
		}

	/* Copy the Work Array Back Into the Code Array */
		for (uint64 i = 0; i < ClosedChains; ++i)
		{
			CodeArray[i] = WorkArray[i];
		}
	}

	cout << "done.\n";

/* *** Step #3: Eliminate Duplicates From List */
	cout << "Eliminate Duplicates ... ";

/* Counter for Unique Polygons (up to reflections) */
	uint64 UniquePolygons = 0;

/* Choic Guarantees First Polygon To Be Considered Unique */
	uint64 PreviousCode = ClosedChains;

/* Collect Unique Polygons */
	for (uint64 i = 0; i < ClosedChains; ++i)
	{
		if (CodeArray[i] != LastCode)
		{
		/* Found a New Type of Polygon */
			LastCode = CodeArray[i];

			WorkArray[UniquePolygons] = LastCode;

			++UniquePolygons;
		}
	}

/* Get Rid of the Code Array */
	delete[] CodeArray;

/* Re-Build Array for Unique Polygon Codes */
	PolyMath *PrimitivePolygon = new PolyMath [UniquePolygons];

/* Copy Results Back Into Polygon Array */
	for (uint64 i = 0; i < UniquePolygons; ++i)
	{
		PrimitivePolygon[i] = PolyMath(WorkArray[i], Length);
	}

/* Memory Clean-Up */
	delete[] WorkArray;

	cout << "done.\n";

/* Step #4:  Examine Symmetry Properties of Chains */
/* (Note:  Possibilities are 1, 2, 3, 6-fold rotational symmetry, perhaps with additional mirror symmetry.) */
	cout << "Examine Symmetry Properties ... ";

/* Set Up Counters For Symmetry Classes */
	uint64 SC1 = 0;
	uint64 SC1m = 0;
	uint64 SC2 = 0;
	uint64 SC2m = 0;
	uint64 SC3 = 0;
	uint64 SC3m = 0;
	uint64 SC6 = 0;
	uint64 SC6m = 0;

/* Analyze Primitive Polygons */
	for (uint64 i = 0; i < UniquePolygons; ++i)
	{
	/* Analyze for Rotational and Mirror Symmetry */
		int RotSym = PrimitivePolygon[i].RotationalSymmetry();
		bool MirrSymm = PrimitivePolygon[i].MirrorSymmetry();

	/* Count Occurrences */
		if (MirrSymm == true)
		{
		/* Symmetry Classes With Mirror Symmetry */
			switch (RotSym)
			{
			case 1:
				++SC1m;
				break;
			case 2:
				++SC2m;
				break;
			case 3:
				++SC3m;
				break;
			case 6:
				++SC6m;
				break;
			default:
				break;
			}
		}
		else
		{
		/* Symmetry Classes Without Mirror Symmetry */
			switch (RotSym)
			{
			case 1:
				++SC1;
				break;
			case 2:
				++SC2;
				break;
			case 3:
				++SC3;
				break;
			case 6:
				++SC6;
				break;
			default:
				break;
			}
		}
	}

	cout << "done.\n\n";

/* Timing Support - End of Calculation */
	FinishTime = clock();
	StartToFinish = Duration(StartTime, FinishTime);

/* Send a Brief Message */
	cout << "(Found " << UniquePolygons << " unique self-avoiding polygon(s) in " << StartToFinish << " seconds) \n\n";

/* For Now, Just Report Results: */

/* Display Final Results */
	cout << "\n *** RESULTS for Chains on 2D Honeycomb Lattice with " << Length << " Segments:\n\n";

/* Number of Non-Overlapping Chains */
	cout << "Number of Non-Overlapping Chains: " << NonOverlaps << "\n\n";

/* Number of Closed-Loop Chains */
	cout << "Number of Closed-Loop Chains: " << ClosedChains << "\n\n";

/* Number of Unique Polygons */
	cout << "Number of Unique Polygons: " << UniquePolygons << " (includes mirror symmetric pairs)\n\n";

/* Specify Polygons By Symmetry Class: */
	cout << "Self-Avoiding Polygon(s) By Symmetry Class: \n\n";

	cout << "Class 1  (trivial symmetry group) ............ " << SC1 << "\n"
		 << "Class 1m (only mirror symmetry) .............. " << SC1m << "\n"
		 << "Class 2  (symmetry under 180° rotations) ..... " << SC2 << "\n"
		 << "Class 2m (180° rotation & mirror symmetry) ... " << SC2m << "\n"
		 << "Class 3  (symmetry under 120° rotations) ..... " << SC3 << "\n"
		 << "Class 3m (120° rotation & mirror symmetry) ... " << SC3m << "\n"
		 << "Class 6  (symmetry under 60° rotations) ...... " << SC6 << "\n"
		 << "Class 6m (60° rotation & mirror symmetry) .... " << SC6m << "\n\n";

/* Clean-Up Memory */
	delete[] ChainArray;
	delete[] PrimitivePolygon;

	return 0;
}