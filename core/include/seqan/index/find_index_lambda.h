// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Enrico Siragusa <enrico.siragusa@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_FIND_INDEX_LAMBDA_H
#define SEQAN_FIND_INDEX_LAMBDA_H

namespace seqan {

// ============================================================================
// Metafunction
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction DefaultFind<Index>
// ----------------------------------------------------------------------------

template <typename THaystack, typename THaystackSpec, typename TPattern>
struct DefaultFind<Index<THaystack, THaystackSpec>, TPattern>
{
    typedef Backtracking<Exact> Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _findImpl(..., Backtracking<Exact>)
// ----------------------------------------------------------------------------

template <typename TState, typename TIndex, typename TNeedle, typename TThreshold, typename TDelegate, typename TSpec>
SEQAN_FUNC_ENABLE_IF(IsSequence<TNeedle>, void)
_findImpl(TState & /* state */,
          TIndex & index,
          TNeedle const & needle,
          TThreshold /* threshold */,
          TDelegate && delegate,
          Backtracking<Exact, TSpec>)
{
    typedef typename Iterator<TIndex, TopDown<> >::Type TIndexIt;

    TIndexIt indexIt(index);

    if (goDown(indexIt, needle))
    {
        delegate(indexIt, TThreshold());
    }
}

// ----------------------------------------------------------------------------
// Function _findImpl(..., Backtracking<Edit/HammingDistance>)
// ----------------------------------------------------------------------------

template <typename TIndexIt, typename TNeedle, typename TNeedleIt, typename TThreshold, typename TDelegate, typename TDistance>
inline void
_findBacktracking(TIndexIt indexIt, TNeedle const & needle, TNeedleIt needleIt, TThreshold errors, TThreshold threshold, TDelegate && delegate, TDistance)
{
    // Exact case.
    if (errors == threshold)
    {
        if (goDown(indexIt, suffix(needle, position(needleIt, needle))))
            delegate(indexIt, errors);
    }
    // Approximate case.
    else if (errors < threshold)
    {
    	//std::cout << suffix(needle, position(needleIt, needle)) << std::endl;
        // Base case.
        if (atEnd(needleIt, needle))
        {
            delegate(indexIt, errors);
        }
        // Recursive case.
        else
        {
            // Insertion.
            if (IsSameType<TDistance, EditDistance>::VALUE)
            {
                _findBacktracking(indexIt, needle, needleIt + 1,
                                  errors + 1, threshold, delegate, TDistance());
            }

            if (goDown(indexIt))
            {
                do
                {
                	// TODO:cpockrandt: oder auch match?
                    // Mismatch.
                    TThreshold delta = !ordEqual(parentEdgeLabel(indexIt), value(needleIt));
                    _findBacktracking(indexIt, needle, needleIt + 1,
                                      errors + delta, threshold, delegate, TDistance());

                    // Deletion.
                    if (IsSameType<TDistance, EditDistance>::VALUE)
                    {
                        _findBacktracking(indexIt, needle, needleIt, errors + 1,
                                          threshold, delegate, TDistance());
                    }
                }
                while (goRight(indexIt));
            }
        }
    }
}

template <typename TLocalHolder, typename TNeedleIter, typename TIndexIter, typename TNeedle, typename TNeedleIt, typename TThreshold, typename TDistance>
inline void
_findBacktrackingSingle(TLocalHolder &, TNeedleIter &, TIndexIter, TNeedle const &, TNeedleIt, TThreshold, TThreshold, TDistance, unsigned long const (&)[LOOKUP_TABLE_SIZE], bool start = true)
{
	std::cerr << "Not possible!" << std::endl;
	return;
}

template <typename TLocalHolder, typename TNeedleIter, typename TText, typename TIndexConfig, typename TTD, typename TNeedle, typename TNeedleIt, typename TThreshold, typename TDistance>
inline void
_findBacktrackingSingle(TLocalHolder & lH,
		  	  	  TNeedleIter & ni,
		  	  	  Iter<Index<TText, BidirectionalFMIndex<TIndexConfig> >, VSTree<TopDown<TTD> > > & indexIt, // TODO: hier war vorher keine referenz!
				  TNeedle const & needle,
				  TNeedleIt needleIt,
                  TThreshold errors,
                  TThreshold threshold,
                  TDistance,
				  unsigned long const (&tbl)[LOOKUP_TABLE_SIZE],
				  bool start = true)
{
    // Exact case.
    if (errors == threshold)
    {
        if (goDown(indexIt, suffix(needle, position(needleIt, needle))))
        	onFindSingleIndex(lH, ni, indexIt);
    }
    // Approximate case.
    else if (errors < threshold)
    {
    	if (start)
    	{
    		if (goDown(indexIt, infixWithLength(needle, position(needleIt, needle), 5)))
    			_findBacktrackingSingle(lH, ni, indexIt, needle, needleIt + 5,
    					errors, threshold, TDistance(), tbl, false);
    		return;
    	}
        // Base case.
        if (atEnd(needleIt, needle))
        	onFindSingleIndex(lH, ni, indexIt);
        // Recursive case.
        else
        {
            // Insertion.
            if (IsSameType<TDistance, EditDistance>::VALUE)
            {
            	_findBacktrackingSingle(lH, ni, indexIt, needle, needleIt + 1,
                                  errors + 1, threshold, TDistance(), tbl, false);
            }

            if (goDown(indexIt))
            {
                do
                {
                    // Mismatch.
                    TThreshold delta = !ordEqual(parentEdgeLabel(indexIt), value(needleIt));
                    _findBacktrackingSingle(lH, ni, indexIt, needle, needleIt + 1,
                                      errors + delta, threshold, TDistance(), tbl, false);

                    // Deletion.
                    if (IsSameType<TDistance, EditDistance>::VALUE)
                    {
                    	_findBacktrackingSingle(lH, ni, indexIt, needle, needleIt, errors + 1,
                                          threshold, TDistance(), tbl, false);
                    }
                }
                while (goRight(indexIt));
            }
        }
    }
}

template <typename TLocalHolder, typename TNeedleIter, typename TIndexIt, typename TNeedle, typename TNeedleIt, typename TThreshold, typename TDistance>
inline void
_findBacktrackingMultiple(TLocalHolder &, TNeedleIter &, TIndexIt, TNeedle const &, TNeedleIt, TThreshold, TThreshold, TDistance,  unsigned long const (&)[LOOKUP_TABLE_SIZE])
{
	std::cerr << "Not possible!" << std::endl;
	return;
}

template <typename TLocalHolder, typename TNeedleIter, typename TIndexIter, typename TNeedle, typename TNeedleIt>
inline void
recursiveSearch(TLocalHolder & lH, TNeedleIter & ni, TIndexIter & indexIt, unsigned int i, unsigned int j,
		unsigned int d, bool goLeft,  unsigned long const (&preComputedIndexPos)[LOOKUP_TABLE_SIZE], TNeedle needle, TNeedleIt _ni)
{
	if (i != 5 || j != 10)
	{
		if (goLeft)
		{
			TNeedle needle2 = value(ni);
			TNeedleIt needleIt2 = begin(needle2, Standard());

			if(!goDown(indexIt.bwdIter, infixWithLength(needle2, position(needleIt2, needle2)+i, 1)/*value(ni)[i]*/))
				return;
		}
		else
		{
			TNeedle needle2 = value(ni+5);
			TNeedleIt needleIt2 = begin(needle2, Standard());

			if(!goDown(indexIt.fwdIter, infixWithLength(needle2, position(needleIt2, needle2)+j-5-1, 1)/*value(ni+5)[j-5-1]*/))
				return;
		}
	}

	if (d == 5)
	{
		onFindSingleIndex(lH, ni+i, indexIt);
		return;
	}

	/*auto vDesc1 = value(indexIt.fwdIter);
	auto vDesc2 = value(indexIt.bwdIter);

	// go left
	if (5 - d == i || 13 <= j)
	{
		recursiveSearch(lH, ni, indexIt, i-1, j, d+1, true, preComputedIndexPos);
	}

	value(indexIt.fwdIter) = vDesc1;
	value(indexIt.bwdIter) = vDesc2;

	// go right
	if (10 + d == j || i <= 2)
	{
		recursiveSearch(lH, ni, indexIt, i, j+1, d+1, false, preComputedIndexPos);
	}*/

	// go left
	if (5 - d == i || 13 <= j)
	{
		if (5 - d == i || 13 <= j)
		{
			auto vDesc1 = value(indexIt.fwdIter);
			auto vDesc2 = value(indexIt.bwdIter);

			recursiveSearch(lH, ni, indexIt, i-1, j, d+1, true, preComputedIndexPos, needle, _ni);

			value(indexIt.fwdIter) = vDesc1;
			value(indexIt.bwdIter) = vDesc2;
		}
		else
			recursiveSearch(lH, ni, indexIt, i-1, j, d+1, true, preComputedIndexPos, needle, _ni);
	}

	// go right
	if (10 + d == j || i <= 2)
	{
		recursiveSearch(lH, ni, indexIt, i, j+1, d+1, false, preComputedIndexPos, needle, _ni);
	}
}

template <typename TIter, typename TNeedle, typename TNeedleIt>
inline bool
_goDownEfficient(TIter & indexIt, TNeedle & needle, TNeedleIt & needleIt,  unsigned long const (&preComputedIndexPos)[LOOKUP_TABLE_SIZE])
{
	unsigned short offset = 5;
	unsigned short effLength = 3;

	if (length(needle)-offset >= effLength)
	{
		unsigned int pos = 4 * (unsigned) ordValue((AminoAcid) needle[offset+0])
		+ 4 * 10 * (unsigned) ordValue((AminoAcid) needle[offset+1])
		+ 4 * 10 * 10 * (unsigned) ordValue((AminoAcid) needle[offset+2]);
		//+ 4 * 10 * 10 * 10 * (unsigned) ordValue((AminoAcid) needle[offset+3]);
		if (preComputedIndexPos[pos+1] > 0)
		{
			indexIt.fwdIter.vDesc.range.i1 = preComputedIndexPos[pos];
			indexIt.fwdIter.vDesc.range.i2 = preComputedIndexPos[pos+1];
			indexIt.fwdIter.vDesc.repLen = effLength;
			indexIt.bwdIter.vDesc.range.i1 = preComputedIndexPos[pos+2];
			indexIt.bwdIter.vDesc.range.i2 = preComputedIndexPos[pos+3];
			indexIt.bwdIter.vDesc.repLen = effLength;
			if (length(needle)-offset > effLength)
			{
				return goDown(indexIt.fwdIter, suffix(needle, position(needleIt, needle) + effLength + offset));
			}
			else
			{
				indexIt.fwdIter.vDesc.lastChar = needle[offset+effLength-1];
				indexIt.bwdIter.vDesc.lastChar = needle[offset+effLength-1];
			}
			return true;
		}
		return false;
	}
	else
		return goDown(indexIt.fwdIter, needle);
}

template <typename TLocalHolder, typename TNeedleIter, typename TText, typename TIndexConfig, typename TTD, typename TNeedle, typename TNeedleIt, typename TThreshold, typename TDistance>
inline void
_findBacktrackingMultiple(TLocalHolder & lH,
		  	  	  TNeedleIter & ni,
		  	  	  Iter<Index<TText, BidirectionalFMIndex<TIndexConfig> >, VSTree<TopDown<TTD> > > & indexIt, // TODO: hier war vorher keine referenz!
				  TNeedle const & needle1,
				  TNeedleIt needleIt1,
                  TThreshold errors,
                  TThreshold threshold,
                  TDistance,
				  unsigned long const (&preComputedIndexPos)[LOOKUP_TABLE_SIZE])
{/*
	// 1-stepping für seed length = 10
	//goDown(indexIt.fwdIter, suffix(needle1, position(needleIt1, needle1)+5)); // TODO: vorher: +5

	//unsigned letter1 = (unsigned) ordValue((AminoAcid) needle1[7]);
	//unsigned letter2 = (unsigned) ordValue((AminoAcid) needle1[8]);
	//unsigned letter3 = (unsigned) ordValue((AminoAcid) needle1[9]);

	//unsigned int pos2 = letter1 + 24 * letter2 + 24 * 24 * letter3;
			//24 * 24 * 24* ordValue(suffix(needle1, position(needleIt1, needle1)+6)[3]);

	if (goDown(indexIt.fwdIter, suffix(needle1, position(needleIt1, needle1)+5)))
			//_goDownEfficient(indexIt, needle1, needleIt1, preComputedIndexPos))
	{
	    //std::cout << "Start: " << indexIt.fwdIter.vDesc.range.i1 << "\t" << indexIt.fwdIter.vDesc.range.i2 << "\t"
	    	//<< indexIt.bwdIter.vDesc.range.i1 << "\t" << indexIt.bwdIter.vDesc.range.i2 << std::endl;


		recursiveSearch(lH, ni, indexIt, 5, 10, 0, false || true, preComputedIndexPos, needle1, needleIt1);
	}
	//goRoot(indexIt);
	return;*/

	// 5-half-overlapping, d=1

	TNeedle needle2 = value(ni+1);

	if (goDown(indexIt.fwdIter, needle2[0]) &&
		goDown(indexIt.fwdIter, needle2[1]) &&
		goDown(indexIt.fwdIter, needle2[2]) &&
		goDown(indexIt.fwdIter, needle2[3])
	)
	{
		// copy ranges, etc.
		auto vDesc1 = value(indexIt.fwdIter);
		auto vDesc2 = value(indexIt.bwdIter);

		if (
			goDown(indexIt.fwdIter, needle2[4]) &&
			goDown(indexIt.fwdIter, needle2[5]) &&
			goDown(indexIt.fwdIter, needle2[6]) &&
			goDown(indexIt.fwdIter, needle2[7])
		)
		{
			onFindSingleIndex(lH, ni+1, indexIt);
		}

		// paste ranges, etc.
		value(indexIt.fwdIter) = vDesc1;
		value(indexIt.bwdIter) = vDesc2;

		if (
			goDown(indexIt.bwdIter, needle1[3]) &&
			goDown(indexIt.bwdIter, needle1[2]) &&
			goDown(indexIt.bwdIter, needle1[1]) &&
			goDown(indexIt.bwdIter, needle1[0])
		)
		{
			onFindSingleIndex(lH, ni, indexIt);
		}

	}
	goRoot(indexIt);

	return;


    // Exact case.
    if (errors == threshold)
    {
        //ModifiedString<TNeedle, ModReverse> const revNeedle(needle);
        //typedef typename Iterator<ModifiedString<TNeedle, ModReverse> const, Standard>::Type  TRevNeedleIt;
        //TRevNeedleIt revNeedleIt = begin(revNeedle, Standard());

        TNeedle needle2 = value(ni+1);
        //TNeedleIt needleIt2 = begin(needle2, Standard());

        /*if(goDown(indexIt.fwdIter, needle))
        {
			onFindSingleIndex(lH, ni, indexIt);
        }

        goRoot(indexIt);

        if(goDown(indexIt.fwdIter, needle2))
        {
			onFindSingleIndex(lH, ni+1, indexIt);
        }

        goRoot(indexIt);*/

    	/*std::cout << "xxx \t" << needle1 << "\t" << needle2
			<< "\t" << infixWithLength(needle2, position(needleIt2, needle2), 4)
			<< "\t" << suffix(needle2, position(needleIt2, needle2)+4)
			<< "\t" << infixWithLength(needle1, position(needleIt1, needle1), 4)
			<< std::endl;*/

    	if (
    		goDown(indexIt.fwdIter, needle2[0]) &&
			goDown(indexIt.fwdIter, needle2[1]) &&
			goDown(indexIt.fwdIter, needle2[2]) &&
			goDown(indexIt.fwdIter, needle2[3])
		)
        //if(goDown(indexIt.fwdIter, needle2))
        //if (goDown(indexIt.fwdIter, infixWithLength(needle2, position(needleIt2, needle2), 4)) &&
    	//		goDown(indexIt.fwdIter, suffix(needle2, position(needleIt2, needle2)+4)))
    	{
			// copy ranges, etc.
    		auto vDesc1 = value(indexIt.fwdIter);
    		auto vDesc2 = value(indexIt.bwdIter);

			if (
				goDown(indexIt.fwdIter, needle2[4]) &&
				goDown(indexIt.fwdIter, needle2[5]) &&
				goDown(indexIt.fwdIter, needle2[6]) &&
				goDown(indexIt.fwdIter, needle2[7])
			)
			{
				onFindSingleIndex(lH, ni+1, indexIt);
			}

			// paste ranges, etc.
			value(indexIt.fwdIter) = vDesc1;
			value(indexIt.bwdIter) = vDesc2;

			if (
				goDown(indexIt.bwdIter, needle1[3]) &&
				goDown(indexIt.bwdIter, needle1[2]) &&
				goDown(indexIt.bwdIter, needle1[1]) &&
				goDown(indexIt.bwdIter, needle1[0])
			)
			{
				onFindSingleIndex(lH, ni, indexIt);
			}

	    	/*
	    	*/
			//if (goDown(indexIt.fwdIter, infixWithLength(needle1, position(needleIt1, needle1), 4)) &&
			//	goDown(indexIt.fwdIter, infixWithLength(needle2, position(needleIt2, needle2), 4)))
	    	//if(goDown(indexIt.fwdIter, needle1))
	    	//{
			//}

    	}
    	goRoot(indexIt);

    	//goRoot(indexIt);

    	/*bool hit1 = false;
    	bool hit2 = false;

		if (goDown(indexIt.fwdIter, infixWithLength(needle2, position(needleIt2, needle2), 4))

		) {
			//typedef typename Iterator<TNeedle>::Type TStrIter;
			//TStrIter strit = end(needle1) - 1;
			//TStrIter stritBegin = begin(needle1);
			unsigned int z = 0;
			//std::cout << "---\t";
			for (; z < 4; ++z)
			//while (!(strit < stritBegin))
			{
				//std::cout << needle1[4-z-1];
				//std::cout << *strit;
				//--strit;
				if (!goDown(indexIt.bwdIter, needle1[4-z-1]))
					break;
			}
			//std::cout << std::endl;

			if (z == 4)
				hit1 = true;
				//onFindSingleIndex(lH, ni, indexIt);
		}

    	goRoot(indexIt);*/

    	/*if (needle1[4] != needle2[0] || needle1[5] != needle2[1] || needle1[6] != needle2[2] || needle1[7] != needle2[3])
    		std::cout
    			<< needle1[0]
    			<< needle1[1]
    			<< needle1[2]
    			<< needle1[3]
    			<< needle1[4]
    			<< needle1[5]
    			<< needle1[6]
    			<< needle1[7] << "\t"
    			<< needle2[0]
    			<< needle2[1]
    			<< needle2[2]
    			<< needle2[3]
    			<< needle2[4]
    			<< needle2[5]
    			<< needle2[6]
    			<< needle2[7] << std::endl;*/

		/*if (goDown(indexIt, suffix(needle2, position(needleIt2, needle2))))
			delegate2(indexIt, errors);

			//goRoot(indexIt);
			//goDown(indexIt.fwdIter, infixWithLength(needle2, position(needleIt2, needle2), 4));

			if (goDown(indexIt, suffix(needle1, position(needleIt1, needle1))))
				delegate1(indexIt, errors);*/

    	            //goRoot(indexIt);
    	    	//}
    	        //if (goDown(indexIt, suffix(needle1, position(needleIt1, needle1))))
    	            //delegate1(indexIt, errors);
    }

        /*if (goDown(indexIt, suffix(needle, position(needleIt, needle))))
	        onFindSingleIndex(lH, ni, indexIt);

        goRoot(indexIt);
    	if (
			goDown(indexIt.fwdIter, infixWithLength(needle2, position(needleIt2, needle2), 4)) &&
			goDown(indexIt.fwdIter, infixWithLength(needle2, position(needleIt2, needle2)+4, 4))
		)
	        onFindSingleIndex(lH, ni+1, indexIt);*/

        //std::cout << needle << "\t" << revNeedle <</* "\t" << infixWithLength(needle, position(needleIt, needle)+4, 4) << "\t"
        		//<< infixWithLength(revNeedle, position(needleIt, needle)+4, 4) <<*/ std::endl;

        /*if (
			goDown(indexIt.fwdIter, needle) &&
			goDown(indexIt.bwdIter, infixWithLength(revNeedle, position(revNeedleIt, revNeedle)+4, 4))
        )
        {


        }
        goRoot(indexIt);
        TNeedle needle2 = value(ni+1);
        TNeedleIt needleIt2 = begin(needle2, Standard());
    	if (
			goDown(indexIt.fwdIter, infixWithLength(needle2, position(needleIt2, needle2), 4)) &&
			goDown(indexIt.fwdIter, infixWithLength(needle2, position(needleIt2, needle2)+4, 4))
		)
	        onFindSingleIndex(lH, ni+1, indexIt);
    }*/
}

/*template <typename TLocalHolder, typename TNeedleIter, typename TState, typename TText, typename TIndexConfig, typename TNeedle,
          typename TThreshold, typename TDistance, typename TSpec>
SEQAN_FUNC_ENABLE_IF(IsSequence<TNeedle>, void)
_findImpl(TLocalHolder &,
		  TNeedleIter &,
		  TState &,
		  Index<TText, BidirectionalFMIndex<TIndexConfig> > & index,
          TNeedle const &,
          TThreshold,
          Backtracking<TDistance, TSpec>)
{
    typedef typename Iterator<Index<TText, BidirectionalFMIndex<TIndexConfig> >, TopDown<> >::Type       TIndexIt;

    TIndexIt indexIt(index);

	goDown(indexIt.bwdIter, "SRVVVNEL");
    std::cout << indexIt.fwdIter.vDesc.range.i1 << "\t" << indexIt.fwdIter.vDesc.range.i2 << "\t"
    	<< indexIt.bwdIter.vDesc.range.i1 << "\t" << indexIt.bwdIter.vDesc.range.i2 << std::endl;
}*/

template <typename TLocalHolder, typename TNeedleIter, typename TState, typename TIndex, typename TNeedle,
          typename TThreshold, typename TDistance, typename TSpec>
SEQAN_FUNC_ENABLE_IF(IsSequence<TNeedle>, void)
_findImplMultiple(TLocalHolder & lH,
		  TNeedleIter & ni,
		  TState & /*indexIt*/,
          TIndex & index,
          TNeedle const & needle,
          TThreshold threshold,
          Backtracking<TDistance, TSpec>,
		  unsigned long const (&preComputedIndexPos)[LOOKUP_TABLE_SIZE])
{
    typedef typename Iterator<TIndex, TopDown<> >::Type       TIndexIt;
    typedef typename Iterator<TNeedle const, Standard>::Type  TNeedleIt;

    TIndexIt indexIt(index);
    TNeedleIt needleIt = begin(needle, Standard());

    _findBacktrackingMultiple(lH, ni, indexIt, needle, needleIt, 0, threshold, TDistance(), preComputedIndexPos);
}

template <typename TLocalHolder, typename TNeedleIter, typename TState, typename TIndex, typename TNeedle,
          typename TThreshold, typename TDistance, typename TSpec>
SEQAN_FUNC_ENABLE_IF(IsSequence<TNeedle>, void)
_findImplSingle(TLocalHolder & lH,
		  TNeedleIter & ni,
		  TState & /*indexIt*/,
          TIndex & index,
          TNeedle const & needle,
          TThreshold threshold,
          Backtracking<TDistance, TSpec>,
		  unsigned long const (&preComputedIndexPos)[LOOKUP_TABLE_SIZE])
{
    typedef typename Iterator<TIndex, TopDown<> >::Type       TIndexIt;
    typedef typename Iterator<TNeedle const, Standard>::Type  TNeedleIt;

    TIndexIt indexIt(index);
    TNeedleIt needleIt = begin(needle, Standard());

    _findBacktrackingSingle(lH, ni, indexIt, needle, needleIt, 0, threshold, TDistance(), preComputedIndexPos, true);
}

template <typename TState, typename TIndex, typename TNeedle,
          typename TThreshold, typename TDelegate, typename TDistance, typename TSpec>
SEQAN_FUNC_ENABLE_IF(IsSequence<TNeedle>, void)
_findImpl(TState & /* indexIt */,
          TIndex & index,
          TNeedle const & needle,
          TThreshold threshold,
          TDelegate && delegate,
          Backtracking<TDistance, TSpec>)
{
    typedef typename Iterator<TIndex, TopDown<> >::Type       TIndexIt;
    typedef typename Iterator<TNeedle const, Standard>::Type  TNeedleIt;

    TIndexIt indexIt(index);
    TNeedleIt needleIt = begin(needle, Standard());
    TThreshold errors = 0;

    _findBacktracking(indexIt, needle, needleIt, errors, threshold, delegate, TDistance());
}

// ----------------------------------------------------------------------------
// Function find(index, index, errors, [](...){}, Backtracking<TDistance>());
// ----------------------------------------------------------------------------

template <typename THaystack, typename THaystackSpec, typename TNeedle, typename TNeedleSpec,
          typename TThreshold, typename TDelegate, typename TDistance, typename TSpec>
inline void
find(Index<THaystack, THaystackSpec> & text,
     Index<TNeedle, TNeedleSpec> & pattern,
     TThreshold threshold,
     TDelegate && delegate,
     Backtracking<TDistance, TSpec>)
{
    typedef typename If<IsSameType<TDistance, Exact>,
                        HammingDistance, TDistance>::Type   TDistance_;
    typedef Backtracking<TDistance_, TSpec>                 TAlgorithm;
    typedef Index<THaystack, THaystackSpec>                 TText;
    typedef Index<TNeedle, TNeedleSpec>                     TPattern;
    typedef Finder_<TText, TPattern, TAlgorithm>            TFinder;

    TFinder finder;

    // This lambda is stored locally because _find() does not accept rvalue delegates.
    std::function<void(TFinder const &)> delegator = [&](TFinder const &)
    {
        delegate(_textIterator(finder), _patternIterator(finder), _getScore(finder));
    };

    _find(finder, text, pattern, threshold, delegator);
}

// ----------------------------------------------------------------------------
// Function find(index, wotdIndex, errors, [](...){}, Backtracking<TDistance>());
// ----------------------------------------------------------------------------

template <typename THaystack, typename THaystackSpec, typename TNeedle, typename TNeedleSpec,
          typename TThreshold, typename TDelegate, typename TDistance, typename TSpec>
inline void
find(Index<THaystack, THaystackSpec> & text,
     Index<TNeedle, IndexWotd<TNeedleSpec> > & pattern,
     TThreshold threshold,
     TDelegate && delegate,
     Backtracking<TDistance, TSpec>)
{
    typedef Index<THaystack, THaystackSpec>                 TText;
    typedef Index<TNeedle, IndexWotd<TNeedleSpec> >         TPattern;
    typedef Backtracking<TDistance, TSpec>                  TAlgorithm;

    typedef Finder<TText, TAlgorithm>                       TFinder_;
    typedef Pattern<TPattern, TAlgorithm>                   TPattern_;

    TFinder_ _finder(text);
    TPattern_ _pattern(pattern, length(back(host(pattern))));

    while (_resume(_finder, _pattern, threshold))
    {
        delegate(_finder.index_iterator, _pattern.index_iterator, _pattern.prefix_aligner.errors);
    }
}

}

#endif  // #ifndef SEQAN_FIND_INDEX_LAMBDA_H
