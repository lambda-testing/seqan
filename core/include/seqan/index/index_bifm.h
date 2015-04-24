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

#ifndef INDEX_BIFM_H_
#define INDEX_BIFM_H_

namespace seqan {

// ==========================================================================
// Forwards
// ==========================================================================

template <typename TSpec = void, typename TConfig = FMIndexConfig<TSpec> >
class BidirectionalFMIndex;

// ==========================================================================
// Classes
// ==========================================================================

// ----------------------------------------------------------------------------
// Class BidirectionalFMIndex
// ----------------------------------------------------------------------------

/*!
 * @class BidirectionalFMIndex
 * @extends Index
 * @headerfile <seqan/index.h>
 * @brief A bidirectional index based on the Burrows-Wheeler transform.
 *
 * @signature template <typename TText[, typename TSpec[, typename TConfig]]>
 *            class Index<TText, FMIndex<TSpec, TConfig> >;
 *
 * @tparam TText   The text type. Types: @link String @endlink, @link StringSet @endlink
 * @tparam TSpec   FM index specialisation, defaults to <tt>void</tt>.
 * @tparam TConfig A config object which determines the data types of the different fibres, defaults to
 *                 <tt>FMIndexConfig&lt;TSpec&gt;</tt>. TBidirectional is only used internally to
 *                 distinguish between a stand-alone FM index and an FM index belonging to a pair of
 *                 FM indices forming a bidirectional FM index.
 *
 * @section Structure
 *
 * The FM index consists of various @link Fibre @endlink of which the most important ones are the compressed
 * suffix array and the LF table, which provides all necessary information for the LF mapping.
 */

template <typename TText, typename TIrgendwas, typename TSpec, typename TSpec2, typename TBidirectional>
class Index<StringSet<TText, TIrgendwas>, BidirectionalFMIndex<TSpec, FMIndexConfig<TSpec2, TBidirectional> > >
{
	typedef ModifiedString<TText, ModReverse>				TRevText;
	typedef Index<StringSet<TRevText, TIrgendwas>, FMIndex<TSpec, FMIndexConfig<TSpec2, FMBidirectional> > >		TRevIndex;
	typedef Index<StringSet<TText, TIrgendwas>	 , FMIndex<TSpec, FMIndexConfig<TSpec2, FMBidirectional> > >			TFwdIndex;

	public:

	StringSet<TRevText, TIrgendwas>	revText;
	TRevIndex		rev;
	TFwdIndex		fwd;

	Index()	{}

	Index(StringSet<TText, TIrgendwas> & text) :
		revText(text),
		rev(revText),
		fwd(text)
	{ }
};

template <typename TText, typename TSpec, typename TSpec2, typename TBidirectional>
class Index<TText, BidirectionalFMIndex<TSpec, FMIndexConfig<TSpec2, TBidirectional> > >
{
	typedef ModifiedString<TText, ModReverse>				TRevText;
	typedef Index<TRevText, FMIndex<TSpec, FMIndexConfig<TSpec2, FMBidirectional> > >		TRevIndex;
	typedef Index<TText, FMIndex<TSpec, FMIndexConfig<TSpec2, FMBidirectional> > >			TFwdIndex;

	public:

	TRevText	revText;
	TRevIndex		rev;
	TFwdIndex		fwd;

	Index()	{}

	Index(TText & text) :
		revText(text),
		rev(revText),
		fwd(text)
	{ }
};

/*template <typename TText, typename TSpec, typename TSpec2, typename TBidirectional>
class Index<TText, BidirectionalFMIndex<TSpec, FMIndexConfig<TSpec2, TBidirectional> > >
{
	typedef ModifiedString<TText, ModReverse>				TRevText;
	typedef Index<TRevText, FMIndex<TSpec, FMIndexConfig<TSpec2, FMBidirectional> > >		TRevIndex;
	typedef Index<TText, FMIndex<TSpec, FMIndexConfig<TSpec2, FMBidirectional> > >			TFwdIndex;

	public:

	TRevText		revText;
	TRevIndex		rev;
	TFwdIndex		fwd;

	Index()	{}

	Index(TText & text) :
		revText(text),
		rev(revText),
		fwd(text)
	{}
};*/

// TODO:cpockrandt: actually we don't need them for FM indices ...
template <typename TText, typename TSpec, typename TConfig>
inline typename Fibre<Index<TText, FMIndex<TSpec, TConfig> >, FibreText>::Type & indexText(Index<TText, BidirectionalFMIndex<TSpec, TConfig> > &index) { return indexText(index.fwd); }
template <typename TText, typename TSpec, typename TConfig>
inline typename Fibre<Index<TText, FMIndex<TSpec, TConfig> > const, FibreText>::Type & indexText(Index<TText, BidirectionalFMIndex<TSpec, TConfig> > const &index) { return indexText(index.fwd); }

template <typename TText, typename TSpec, typename TSpec2, typename TBidirectional>
SEQAN_HOST_DEVICE inline typename Fibre<Index<TText, FMIndex<TSpec, FMIndexConfig<TSpec2, FMBidirectional> > >, FibreSA>::Type & indexSA(Index<TText, BidirectionalFMIndex<TSpec, FMIndexConfig<TSpec2, TBidirectional> > > &index) { return indexSA(index.fwd); }
template <typename TText, typename TSpec, typename TSpec2, typename TBidirectional>
SEQAN_HOST_DEVICE inline typename Fibre<Index<TText, FMIndex<TSpec, FMIndexConfig<TSpec2, FMBidirectional> > > const, FibreSA>::Type & indexSA(Index<TText, BidirectionalFMIndex<TSpec, FMIndexConfig<TSpec2, TBidirectional> > > const &index) { return indexSA(index.fwd); }


// This function can be used to open a previously saved index.
template <typename TText, typename TSpec, typename TConfig>
inline bool open(Index<TText, BidirectionalFMIndex<TSpec, TConfig> > & index, const char * fileName)
{
	String<char> name;

	name = fileName;    append(name, ".fwd");
    bool fwdIndex = open(index.fwd, toCString(name), DefaultOpenMode<Index<TText, FMIndex<TSpec, TConfig> > >::VALUE);
    if (fwdIndex)
    {
    	name = fileName;    append(name, ".bwd");
    	return open(index.rev, toCString(name), DefaultOpenMode<Index<TText, FMIndex<TSpec, TConfig> > >::VALUE);
    }
    return false;
}

// This function can be used to save an index on disk.
template <typename TText, typename TSpec, typename TConfig>
inline bool save(Index<TText, BidirectionalFMIndex<TSpec, TConfig> > const & index, const char * fileName)
{
	String<char> name;

	name = fileName;    append(name, ".fwd");
	bool fwdIndex = save(index.fwd, toCString(name), DefaultOpenMode<Index<TText, FMIndex<TSpec, TConfig> > >::VALUE);
	if (fwdIndex)
	{
		name = fileName;    append(name, ".bwd");
		return save(index.rev, toCString(name), DefaultOpenMode<Index<TText, FMIndex<TSpec, TConfig> > >::VALUE);
	}
    return false;
}

}
#endif /* INDEX_BIFM_H_ */
