// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2015, Knut Reinert, FU Berlin
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
// Author: Christopher Pockrandt <christopher.pockrandt@fu-berlin.de>
// ==========================================================================

//SEQAN_NO_DDDOC:do not generate documentation for this file

#ifndef INDEX_BIFM_STREE_H_
#define INDEX_BIFM_STREE_H_

namespace seqan {

// ==========================================================================
// Classes
// ==========================================================================

// ----------------------------------------------------------------------------
// Class BidirectionalFMIndex-Iter
// ----------------------------------------------------------------------------
//TODO:cpockrandt:documentation
template <typename TText, class TOccSpec, typename TBidirectional, typename TSpec>
class Iter<Index<TText, BidirectionalFMIndex<TOccSpec, FMIndexConfig<TOccSpec, TBidirectional> > >, VSTree<TopDown<ParentLinks<TSpec> > > >
{
public:
	typedef Index<TText, BidirectionalFMIndex<TOccSpec, FMIndexConfig<TOccSpec, TBidirectional> > >	TBiIndex;

	typedef ModifiedString<TText, ModReverse>																	TRevText;

	typedef Index<TText, FMIndex<TOccSpec, FMIndexConfig<TOccSpec, FMBidirectional> > >				TFwdIndex;
	typedef Index<TRevText, FMIndex<TOccSpec, FMIndexConfig<TOccSpec, FMBidirectional> > >			TRevIndex;

	typedef Iter<TFwdIndex, VSTree<TopDown<ParentLinks<TSpec> > > >												TFwdIndexIter;
	typedef Iter<TRevIndex, VSTree<TopDown<ParentLinks<TSpec> > > >												TRevIndexIter;

	TFwdIndexIter	fwdIter;
	TRevIndexIter	bwdIter;

//____________________________________________________________________________

	Iter():
		fwdIter(),
		bwdIter() {}

	Iter(TBiIndex &_index):
		fwdIter(*(&_index.fwd)),
		bwdIter(*(&_index.rev))
	{
		fwdIter.setRevIter(bwdIter);
		bwdIter.setRevIter(fwdIter);
	}

	Iter(TBiIndex &_index, MinimalCtor):
		fwdIter(*(&_index.fwd), MinimalCtor()),
		bwdIter(*(&_index.rev), MinimalCtor())
	{
		fwdIter.setRevIter(bwdIter);
		bwdIter.setRevIter(fwdIter);
	}

	template <typename TSpec2>
	Iter(Iter<TBiIndex, VSTree<TopDown<TSpec2> > > const &_origin):
		fwdIter(*(&_origin.fwdIter)),
		bwdIter(*(&_origin.bwdIter))
	{
		fwdIter.setRevIter(bwdIter);
		bwdIter.setRevIter(fwdIter);
	}

//____________________________________________________________________________

	template <typename TSpec2>
	inline Iter const &
	operator = (Iter<TBiIndex, VSTree<TopDown<TSpec2> > > const &_origin)
	{
		fwdIter = *(&_origin.fwdIter);
		bwdIter = *(&_origin.bwdIter);
		return *this;
	}
};




//----------------------------------------------------------------------
//----------------------------------------------------------------------
//----------------------------------------------------------------------
//----------------------------------------------------------------------




template <typename TText, class TOccSpec, typename TSpec>
class Iter<Index<ModifiedString<TText, ModReverse>, FMIndex<TOccSpec, FMIndexConfig<TOccSpec, FMBidirectional> > >, VSTree<TopDown<ParentLinks<TSpec> > > >
{
public:

	typedef ModifiedString<TText, ModReverse>															TRevText;
	typedef Index<TText, FMIndex<TOccSpec, FMIndexConfig<TOccSpec, FMBidirectional> > >		TRevIndex;
	typedef Index<TRevText, FMIndex<TOccSpec, FMIndexConfig<TOccSpec, FMBidirectional> > >	TFwdIndex;

	typedef typename VertexDescriptor<TFwdIndex>::Type	TVertexDesc;

	typedef Iter										iterator;
	typedef Iter<TFwdIndex, VSTree<TopDown<ParentLinks<TSpec> > > >	TFwdIndexIter;
	typedef Iter<TRevIndex, VSTree<TopDown<ParentLinks<TSpec> > > >	TRevIndexIter;

	typedef	typename HistoryStackEntry_<Iter>::Type		TStackEntry;
	typedef String<TStackEntry, Block<> >				TStack;

	TFwdIndex const	*index;		// container of all necessary tables of the forward iterator
	TRevIndexIter 	*revIter;	// container of all necessary tables of the backward iterator
	TVertexDesc		vDesc;		// current interval in suffix array and
								// right border of parent interval (needed in goRight)

	// pseudo history stack (to go up at most one node)
	TVertexDesc		_parentDesc;

	TStack			history;	// contains all previously visited intervals (allows to go up)

//____________________________________________________________________________

	Iter() : index() {}

	Iter(TFwdIndex &_index):
		index(&_index)
	{
		_indexRequireTopDownIteration(_index);
		goRoot(*this);
	}

	Iter(TFwdIndex &_index, MinimalCtor):
		index(&_index),
		vDesc(MinimalCtor()),
		_parentDesc(MinimalCtor()) {}

	// NOTE(esiragusa): _parentDesc is unitialized
	Iter(TFwdIndex &_index, TVertexDesc const &_vDesc):
		index(&_index),
		vDesc(_vDesc)
	{
		_indexRequireTopDownIteration(_index);
	}

	template <typename TSpec2>
	Iter(Iter<TFwdIndex, VSTree<TopDown<TSpec2> > > const &_origin):
		index(&container(_origin)),
		vDesc(value(_origin)),
		_parentDesc(nodeUp(_origin))/*,
		history(_origin.histor)*/ {}

//____________________________________________________________________________

	template <typename TSpec2>
	inline Iter const &
	operator = (Iter<TFwdIndex, VSTree<TopDown<TSpec2> > > const &_origin)
	{
		index = &container(_origin);
		vDesc = value(_origin);
		_parentDesc = nodeUp(_origin);
		//history = _origin.history;
		return *this;
	}

	void setRevIter(TRevIndexIter &_revIter)
	{
		revIter = &_revIter;
	}
};

template <typename TText, class TOccSpec, typename TSpec>
class Iter<Index<TText, FMIndex<TOccSpec, FMIndexConfig<TOccSpec, FMBidirectional> > >, VSTree<TopDown<ParentLinks<TSpec> > > >
{
public:

	typedef Iter	iterator;

	typedef Index<TText, FMIndex<TOccSpec, FMIndexConfig<TOccSpec, FMBidirectional> > >		TFwdIndex;
	typedef ModifiedString<TText, ModReverse>															TRevText;
	typedef Index<TRevText, FMIndex<TOccSpec, FMIndexConfig<TOccSpec, FMBidirectional> > >	TRevIndex;

	typedef Iter<TFwdIndex, VSTree< TopDown<ParentLinks<TSpec> > > >	TFwdIndexIter;
	typedef Iter<TRevIndex, VSTree< TopDown<ParentLinks<TSpec> > > >	TRevIndexIter;

	typedef typename VertexDescriptor<TFwdIndex>::Type	TVertexDesc;

	typedef	typename HistoryStackEntry_<Iter>::Type		TStackEntry;
	typedef String<TStackEntry, Block<> >				TStack;

	TFwdIndex const	*index;		// container of all necessary tables
	TRevIndexIter 	*revIter;	// container of all necessary tables
	TVertexDesc		vDesc;		// current interval in suffix array and
								// right border of parent interval (needed in goRight)

	// pseudo history stack (to go up at most one node)
	TVertexDesc		_parentDesc;

	TStack			history;	// contains all previously visited intervals (allows to go up)

//____________________________________________________________________________

	Iter() : index() {}

	Iter(TFwdIndex &_index):
		index(&_index)
	{
		_indexRequireTopDownIteration(_index);
		goRoot(*this);
	}

	Iter(TFwdIndex &_index, MinimalCtor):
		index(&_index),
		vDesc(MinimalCtor()),
		_parentDesc(MinimalCtor()) {}

	// NOTE(esiragusa): _parentDesc is unitialized
	Iter(TFwdIndex &_index, TVertexDesc const &_vDesc):
		index(&_index),
		vDesc(_vDesc)
	{
		_indexRequireTopDownIteration(_index);
	}

	template <typename TSpec2>
	Iter(Iter<TFwdIndex, VSTree<TopDown<TSpec2> > > const &_origin):
		index(&container(_origin)),
		vDesc(value(_origin)),
		_parentDesc(nodeUp(_origin)),
		history(_origin.history) {}

//____________________________________________________________________________

	template <typename TSpec2>
	inline Iter const &
	operator = (Iter<TFwdIndex, VSTree<TopDown<TSpec2> > > const &_origin)
	{
		index = &container(_origin);
		vDesc = value(_origin);
		_parentDesc = nodeUp(_origin);
		history = _origin.history;
		return *this;
	}

	void setRevIter(TRevIndexIter &_revIter)
	{
		revIter = &_revIter;
	}
};

/*template <typename TText, class TOccSpec, class TLengthSum, typename TSpec>
class Iter<Index<TText, FMIndex<TOccSpec, FMIndexConfig<TOccSpec, TLengthSum, FMBidirectional> > >, VSTree<TopDown<ParentLinks<TSpec> > > >:
    public Iter<Index<TText, FMIndex<TOccSpec, FMIndexConfig<TOccSpec, TLengthSum, FMBidirectional> > >, VSTree<TopDown<> > >
{
public:
	typedef Index<TText, FMIndex<TOccSpec, FMIndexConfig<TOccSpec, TLengthSum, FMBidirectional> > >	TIndex;
    typedef Iter<TIndex, VSTree<TopDown<> > >	TBase;
    typedef typename HistoryStack_<Iter>::Type	TStack;
    typedef Iter								iterator;

    TStack            history;    // contains all previously visited intervals (allows to go up)

//____________________________________________________________________________

    SEQAN_HOST_DEVICE
    Iter() :
        TBase() {}

    SEQAN_HOST_DEVICE
    Iter(TIndex &_index):
        TBase(_index) {}

    SEQAN_HOST_DEVICE
    Iter(TIndex &_index, MinimalCtor):
        TBase(_index, MinimalCtor()) {}

    SEQAN_HOST_DEVICE
    Iter(Iter const &_origin):
        TBase((TBase const &)_origin),
        history(_origin.history) {}

//____________________________________________________________________________

    SEQAN_HOST_DEVICE inline
    Iter const &
    operator = (Iter const &_origin)
    {
        *(TBase*)(this) = _origin;
        history = _origin.history;
        return *this;
    }
};

template <typename TText, class TOccSpec, class TLengthSum, typename TSpec>
class Iter<Index<ModifiedString<TText, ModReverse>, FMIndex<TOccSpec, FMIndexConfig<TOccSpec, TLengthSum, FMBidirectional> > >, VSTree<TopDown<ParentLinks<TSpec> > > >:
    public Iter<Index<TText, FMIndex<TOccSpec, FMIndexConfig<TOccSpec, TLengthSum, FMBidirectional> > >, VSTree<TopDown<> > >
{
public:
	typedef Index<ModifiedString<TText, ModReverse>, FMIndex<TOccSpec, FMIndexConfig<TOccSpec, TLengthSum, FMBidirectional> > >	TIndex;
    typedef Iter<TIndex, VSTree<TopDown<> > >	TBase;
    typedef typename HistoryStack_<Iter>::Type	TStack;
    typedef Iter								iterator;

    TStack            history;    // contains all previously visited intervals (allows to go up)

//____________________________________________________________________________

    SEQAN_HOST_DEVICE
    Iter() :
        TBase() {}

    SEQAN_HOST_DEVICE
    Iter(TIndex &_index):
        TBase(_index) {}

    SEQAN_HOST_DEVICE
    Iter(TIndex &_index, MinimalCtor):
        TBase(_index, MinimalCtor()) {}

    SEQAN_HOST_DEVICE
    Iter(Iter const &_origin):
        TBase((TBase const &)_origin),
        history(_origin.history) {}

//____________________________________________________________________________

    SEQAN_HOST_DEVICE inline
    Iter const &
    operator = (Iter const &_origin)
    {
        *(TBase*)(this) = _origin;
        history = _origin.history;
        return *this;
    }
};*/







// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _getNodeByChar()                                         [Iterator]
// ----------------------------------------------------------------------------

// TODO:christopher: simplify, documentation, commands
template <typename TText, typename TOccSpec, class TSpec, typename TChar>
inline bool _getNodeByChar(Iter<Index<TText, FMIndex<TOccSpec, FMIndexConfig<TOccSpec, FMBidirectional> > >, VSTree<TopDown<TSpec> > > const & it,
                           typename VertexDescriptor<Index<TText, FMIndex<TOccSpec, FMIndexConfig<TOccSpec, FMBidirectional> > > >::Type const & vDesc,
                           Pair<typename Size<Index<TText, FMIndex<TOccSpec, FMIndexConfig<TOccSpec, FMBidirectional> > > >::Type> & _range,
                           TChar c)
{
    typedef Index<TText, FMIndex<TOccSpec, FMIndexConfig<TOccSpec, FMBidirectional> > >        TIndex;
    typedef typename Value<TIndex>::Type                        TAlphabet;

	typedef typename Fibre<TIndex, FibreLF>::Type               TLF;

	Pair<typename Size<Index<TText, FMIndex<TOccSpec, FMIndexConfig<TOccSpec, FMBidirectional> > > >::Type> tmpRange;

	TIndex const & index = container(it);
	TLF const & lf = indexLF(index);

	unsigned int sum = 0;
	int alpSize = ValueSize<TAlphabet>::VALUE;

	// TODO: int nehmen statt TAlphabet?
	for (TAlphabet _d = 0; _d < std::min(alpSize, (int) c); ++_d)
	{
		unsigned int i1, i2;
		TChar d = TAlphabet(_d);

		tmpRange = range(index, vDesc);
		i1 = lf(tmpRange.i1, d);
		i2 = lf(tmpRange.i2, d);

		if (i1 < i2)
		{
			sum += i2 - i1;
		}
	}

	// TODO:chris
	unsigned int sentinelPos = lf.sentinels;//index.lfTable.occTable.sentinelPosition;

	if (/*vDesc.repLen > 0 && */vDesc.range.i1 <= sentinelPos && sentinelPos < vDesc.range.i2)
		++sum;

	bool isRoot = _isRoot(vDesc);
	_range = range(index, vDesc);
	_range.i1 = lf(_range.i1, c);
	_range.i2 = lf(_range.i2, c);

    if (_range.i1 < _range.i2)
    {
    	// historyPush nicht änderbar, wg. index_fm_stree.h Z. 300 ff. und eigener _getNodeByChar-Impl.! Gibt sonst Probleme mit der Reihenfolge
        _historyPush(*it.revIter); // "it" itself is already pushed in the wrapping method

		if (isRoot)
		{
			value(*it.revIter).range.i1 = _range.i1;
			value(*it.revIter).range.i2 = _range.i2;
		}
		else
		{
			value(*it.revIter).range.i1 += sum;
			value(*it.revIter).range.i2 = value(*it.revIter).range.i1 + (_range.i2 - _range.i1);
		}
		value(*it.revIter).lastChar = c;
		value(*it.revIter).repLen++;
        return true;
    }

    return false;
}

// ----------------------------------------------------------------------------
// Function leftExtend()                                         [Iterator]
// ----------------------------------------------------------------------------
// TODO:christopher: generalize for sequence objects, documentation
template <typename TText, class TIndexSpec, class TOccSpec, typename TSpec, typename TChar>
inline
bool leftExtend(Iter<Index<TText, BidirectionalFMIndex<TIndexSpec, TOccSpec> >, VSTree<TopDown<TSpec> > > & it, TChar c)
{
	return goDown(it.fwdIter, c);
}

// ----------------------------------------------------------------------------
// Function leftExtend()                                         [Iterator]
// ----------------------------------------------------------------------------
// TODO:christopher: generalize for sequence objects, documentation
/*template <typename TText, class TIndexSpec, class TOccSpec, typename TSpec, typename TChar>
inline
bool leftExtend(Iter<Index<TText, BidirectionalFMIndex<TIndexSpec, TOccSpec> >, VSTree<TopDown<ParentLinks<TSpec> > > > & it, TText str)
{
	//return goDown(it.fwdIter, str, length(str));
	for (unsigned int i = 0; i < length(str); ++i) {
		if(!goDown(it.fwdIter, str[i]))
			return false;
	}
	return true;
}*/

// ----------------------------------------------------------------------------
// Function rightExtend()                                         [Iterator]
// ----------------------------------------------------------------------------
// TODO:christopher: generalize for sequence objects, documentation
template <typename TText, class TIndexSpec, class TOccSpec, typename TSpec, typename TChar>
inline
bool rightExtend(Iter<Index<TText, BidirectionalFMIndex<TIndexSpec, TOccSpec> >, VSTree<TopDown<TSpec> > > & it, TChar c)
{
	return goDown(it.bwdIter, c);
}

template <typename TText, typename TOccSpec, typename TIndexSpec, typename TSpec>
void goRoot(Iter<Index<TText, BidirectionalFMIndex<TOccSpec, TIndexSpec> >, VSTree<TopDown<TSpec> > > & it)
{
    goRoot(it.fwdIter);
    goRoot(it.bwdIter);
}

/*template <typename TText, typename TOccSpec, typename TIndexSpec, typename TSpec>
bool goUp(Iter<Index<TText, BidirectionalFMIndex<TOccSpec, TIndexSpec> >, VSTree<TopDown<TSpec> > > & it)
{
	//if (isRoot(it)) return false;



	goUp(it.fwdIter);
	return false;



	//goUp(it.bwdIter);// && goUp(it.bwdIter);// && goUp(*it.revIter);
	//return false;
	//value(it).lastChar = back(it.history).lastChar;
	//_historyPop(it);

	//value(*it.revIter).lastChar = back((*it.revIter).history).lastChar;
	//_historyPop((*it.revIter));
	//return true;
}*/

template <typename TText, typename TOccSpec, typename TIndexSpec, typename TSpec>
bool goUp(Iter<Index<TText, BidirectionalFMIndex<TOccSpec, TIndexSpec> >, VSTree<TopDown<ParentLinks<TSpec> > > > & it)
{
	//if (isRoot(it)) return false;



	//goUp(it.fwdIter);
	return true;



	//goUp(it.bwdIter);// && goUp(it.bwdIter);// && goUp(*it.revIter);
	//return false;
	//value(it).lastChar = back(it.history).lastChar;
	//_historyPop(it);

	//value(*it.revIter).lastChar = back((*it.revIter).history).lastChar;
	//_historyPop((*it.revIter));
	//return true;
}

template < typename TText, typename TOccSpec, typename TIndexSpec, class TSpec >
SEQAN_HOST_DEVICE inline typename Infix<typename Fibre<Index<TText, BidirectionalFMIndex<TOccSpec, TIndexSpec> >, FibreSA>::Type const >::Type
getOccurrences(Iter<Index<TText, BidirectionalFMIndex<TOccSpec, TIndexSpec> >, VSTree<TSpec> > const &it)
{
	return getOccurrences(it.fwdIter);
}







/**************************************************************************
 * HISTORY STUFF *
 */





	// TODO: pointer, etc. checken
	/*template < typename TText, class TOccSpec, class TIndexSpec >
	inline void
	_historyClear(Iter<Index<TText, FMIndex<TOccSpec, TIndexSpec> >, VSTree< TopDown<BiDirectional> > > &it)
	{
		it._parentDesc = value(it);
		// when first iterator is initialized, this method is called but revIter is not initialized yet!
		//(*it.revIter)._parentDesc = value(*it.revIter);
	}

	template < typename TText, class TOccSpec, class TIndexSpec >
	inline void
	_historyClear(Iter<Index<TText, FMIndex<TOccSpec, TIndexSpec> >, VSTree< TopDown< ParentLinks<BiDirectional> > > > &it)
	{
		clear(it.history);
		//clear((*it.revIter).history);
	}*/



	/*template <typename TText, typename TOccSpec, typename TIndexSpec>
	inline void _historyPush(Iter<Index<TText, FMIndex<TOccSpec, TIndexSpec > >, VSTree<TopDown<BiDirectional> > > & it)
	{
	    it._parentDesc = value(it);
	    //(*it.revIter)._parentDesc = value(*it.revIter);
	}*/

	/*template < typename TText, typename TOccSpec, typename TIndexSpec>
	inline void
	_historyPush(Iter<Index<TText, FMIndex<TOccSpec, TIndexSpec > >, VSTree<TopDown<ParentLinks<BiDirectional> > > > & it)
	{
	    typedef Iter<Index<TText, FMIndex<TOccSpec, TIndexSpec > >, VSTree<TopDown<ParentLinks<BiDirectional> > > > TIter;

	    typename HistoryStackEntry_<TIter>::Type h;
	    //typename HistoryStackEntry_<TIter>::Type h2;

	    h.range = value(it).range;
	    h.lastChar = value(it).lastChar;
	    //h2.range = value(*it.revIter).range;
	    //h2.lastChar = value(*it.revIter).lastChar;

	    appendValue(it.history, h);
	    //appendValue((*it.revIter).history, h2);
	}*/




	/*template <typename TText, typename TOccSpec, class TLengthSum, typename TSpec>
	bool _goUp(Iter<Index<TText, FMIndex<TOccSpec, FMIndexConfig<TOccSpec, TLengthSum, FMBidirectional> > >, VSTree<TopDown<ParentLinks<TSpec> > > > & it)
	{
		//if (isRoot(it)) return false;

		return _goUp(it) && _goUp(*it.revIter);

		//value(it).lastChar = back(it.history).lastChar;
		//_historyPop(it);

		//value(*it.revIter).lastChar = back((*it.revIter).history).lastChar;
		//_historyPop((*it.revIter));
		//return true;
	}*/

}

#endif /* INDEX_BIFM_STREE_H_ */
