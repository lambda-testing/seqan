// ==========================================================================
//                                  lambda
// ==========================================================================
// Copyright (c) 2013, Hannes Hauswedell, FU Berlin
// All rights reserved.
//
// This file is part of Lambda.
//
// Lambda is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Lambda is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Lambda.  If not, see <http://www.gnu.org/licenses/>.*/
// ==========================================================================
// Author: Hannes Hauswedell <hauswedell@mi.fu-berlin.de>
// ==========================================================================
// store.h: contains types and definitions for storing sequences and indices
// ==========================================================================

#ifndef SEQAN_LAMBDA_MISC_H_
#define SEQAN_LAMBDA_MISC_H_

#include <type_traits>
#include <forward_list>

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include <seqan/index.h>

#include <seqan/align.h>
#include <seqan/blast.h>
// #include <seqan/reduced_aminoacid.h>

#include "options.hpp"

using namespace seqan;

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Metafunctions
// ============================================================================

// makes partial function specialization convenient
template <bool condition>
using MyEnableIf = typename std::enable_if<condition, int>::type;

// ============================================================================
// Functions for translation and retranslation
// ============================================================================

// template <typename TString, typename TSpec>
// inline void
// assign(StringSet<ModifiedString<TString, TSpec>, Owner<ConcatDirect<>> > & target,
//        StringSet<TString, Owner<ConcatDirect<>> > & source)
// {
//     target.limits = source.limits;
//     target.concat._host = &source.concat;
// }


template <typename T>
inline uint64_t
length(std::deque<T> const & list)
{
    return list.size();
}
template <typename T>
inline uint64_t
length(std::forward_list<T> const & list)
{
    return std::distance(list.begin(), list.end());
}

template <typename TAlph>
inline std::basic_ostream<char> &
operator<<(std::basic_ostream<char> & out,
           const Iter<const String<SimpleType<unsigned char,TAlph>,
                                    seqan::Packed<> >,
                      seqan::Packed<> > it)
{
    out << *it;
    return out;
}


template <typename T1, typename T2>
inline uint64_t
quickHamming(T1 const & s1, T2 const & s2)
{
    SEQAN_ASSERT_EQ(length(s1), length(s2));

    uint64_t ret = 0;

    for (uint64_t i = 0; i < length(s1); ++i)
        if (s1[i] != s2[i])
            ++ret;

    return ret;
}

template <typename TPos>
inline bool
inRange(TPos const i, TPos const beg, TPos const end)
{
    return ((i >= beg) && (i < end));
}

inline int64_t
intervalOverlap(uint64_t const s1, uint64_t const e1,
                uint64_t const s2, uint64_t const e2)
{
    return std::min(e1, e2) - std::max(s1, s2);
}

inline void
printProgressBar(uint64_t & lastPercent, uint64_t curPerc)
{
    //round down to even
    curPerc = curPerc & ~1;
//     #pragma omp critical(stdout)
    if ((lastPercent != curPerc) && (curPerc <= 100))
    {
        for (uint64_t i = lastPercent + 2; i <= curPerc; i+=2)
        {
            if (i == 100)
                std::cout << "|\n" << std::flush;
            else if (i % 10 == 0)
                std::cout << "*" << std::flush;
            else
                std::cout << "·" << std::flush;
        }
        lastPercent = curPerc;
    }
}


template <typename TSequence, typename TAlignSpec,
          typename TScoreValue, typename TScoreSpec, typename TAlignContext>
inline TScoreValue
localAlignment2(Align<TSequence, TAlignSpec> & align,
                Score<TScoreValue, TScoreSpec> const & scoringScheme,
                int lowerDiag,
                int upperDiag,
                TAlignContext & alignContext)
{
//     typedef Align<TSequence, TAlignSpec> TAlign;
//     typedef typename Size<TAlign>::Type TSize;
//     typedef typename Position<TAlign>::Type TPosition;
//     typedef TraceSegment_<TPosition, TSize> TTraceSegment;

    SEQAN_ASSERT_EQ(length(rows(align)), 2u);

    clear(alignContext.traceSegment);

    typedef FreeEndGaps_<True, True, True, True> TFreeEndGaps;
    typedef AlignConfig2<LocalAlignment_<>,
                         DPBand,
                         TFreeEndGaps,
                         TracebackOn<TracebackConfig_<CompleteTrace,
                                                      GapsLeft> > > TAlignConfig;

    TScoreValue score;
    DPScoutState_<Default> scoutState;
    score = _setUpAndRunAlignment(alignContext.dpContext,
                                  alignContext.traceSegment,
                                  scoutState,
                                  source(row(align, 0)),
                                  source(row(align, 1)),
                                  scoringScheme,
                                  TAlignConfig(lowerDiag, upperDiag));

    _adaptTraceSegmentsTo(row(align, 0), row(align, 1), alignContext.traceSegment);
    return score;
}


template <typename TText, typename TSpec>
struct ComparisonCounter;

// no counting
template <typename TText>
struct ComparisonCounter<TText, Nothing>
{
    uint64_t _comparisons = 0;
    uint64_t _expectedComparisons = 0;
    uint64_t _lastPercent = 0;
    ComparisonCounter(TText const &,
                      uint64_t expectedComparisons = 0u)
    {
        (void)expectedComparisons;
    }

    // may be constexpr in c++14
    inline void inc() const
    {}
};

// every thread counts
#ifdef _OPENMP
template <typename TText>
struct ComparisonCounter<TText, std::false_type>
#else
template <typename TText, typename TSpec>
struct ComparisonCounter
#endif
{
    uint64_t _comparisons = 0;
    uint64_t _expectedComparisons = 0;
//     uint64_t _twoPercent = 0;
    uint64_t _lastPercent = 0;
    uint64_t _checkEveryNHits = 1;

    ComparisonCounter(TText const & text,
                      uint64_t expectedComparisons = 0u)
    {
        if (expectedComparisons == 0)
        {
            uint64_t l = length(concat(text));
            _expectedComparisons = 1.2 * double(l) * std::log(l) / std::log(2);
        } else
            _expectedComparisons = expectedComparisons;

//         _twoPercent = _expectedComparisons / 50;
        _comparisons = 0;
        _lastPercent = 0;
        while ((_checkEveryNHits << 1) < (_expectedComparisons / 100))
            _checkEveryNHits <<= 1;
    }

    inline void inc()
    {
        uint64_t comp = ++_comparisons;
        // it is not important that the henceforth _comparisons be actually
        // the same value (might not be due to SMP)

        // progress reporting
        if (comp & _checkEveryNHits)
        {
            uint64_t curPerc = comp * 100 / _expectedComparisons;
            if (curPerc < 100)
                printProgressBar(_lastPercent, curPerc);
        }
    }
};

// only one thread counts
#ifdef _OPENMP
template <typename TText>
struct ComparisonCounter<TText, std::true_type>
{
    uint64_t _comparisons = 0;
    uint64_t _expectedComparisons = 0;
//     uint64_t _twoPercent = 0;
    uint64_t _lastPercent = 0;
    uint64_t _checkEveryNHits = 1;

    ComparisonCounter(TText const & text,
                      uint64_t expectedComparisons = 0u)
    {
        if (expectedComparisons == 0)
        {
            uint64_t l = length(concat(text));
            _expectedComparisons = 1.2 * double(l) * std::log(l) / std::log(2) /
                                   omp_get_max_threads();
        } else
            _expectedComparisons = expectedComparisons;

//         _twoPercent = _expectedComparisons / 50;
//         _comparisons = 0;
        while ((_checkEveryNHits << 1) < (_expectedComparisons / 100))
            _checkEveryNHits <<= 1;
    }

    inline void inc()
    {
        if (omp_get_thread_num() == 0) // only one thread counts
        {
            uint64_t comp = ++_comparisons;

            // progress reporting
            if (comp & _checkEveryNHits)
            {
                uint64_t curPerc = comp * 100 / _expectedComparisons;
                if (curPerc < 100)
                    printProgressBar(_lastPercent, curPerc);
            }
        }
    }
};
#endif


template <typename TSA,
          typename TString,
          typename TSSetSpec,
          typename TAlgo,
          typename TLambda>
inline void
createSuffixArray(TSA & SA,
                  StringSet<TString, TSSetSpec> const & s,
                  TAlgo const &,
                  TLambda const &)
{
    return createSuffixArray(SA, s, TAlgo());
}

template <typename TText, typename TSpec, typename TConfig, typename TLambda>
inline bool indexCreate(Index<TText, FMIndex<TSpec, TConfig> > & index,
                        TText const & text,
                        FibreSALF const &,
                        TLambda const & progressCallback)
{
	typedef Index<TText, FMIndex<TSpec, TConfig> >      TIndex;
    typedef typename Fibre<TIndex, FibreTempSA>::Type   TTempSA;
    typedef typename DefaultIndexCreator<TIndex, FibreSA>::Type  TAlgo;
    typedef typename Size<TIndex>::Type                 TSize;

    indexText(index) = text;

    if (empty(text))
        return false;

    TTempSA tempSA;

    // Create the full SA.
    resize(tempSA, lengthSum(text), Exact());
    createSuffixArray(tempSA, text, TAlgo(), progressCallback);

    // Create the LF table.
    createLF(indexLF(index), text, tempSA);

    // Set the FMIndex LF as the CompressedSA LF.
    setFibre(indexSA(index), indexLF(index), FibreLF());

    // Create the compressed SA.
    TSize numSentinel = countSequences(text);
    createCompressedSa(indexSA(index), tempSA, numSentinel);

    return true;
}

//TODO:cpockrandt: wenn was nicht funzt, hier nachsehen!
template <typename TText, typename TIrgendwas, typename TSpec, typename TConfig, typename TLambda>
inline bool indexCreate(Index<StringSet<TText, TIrgendwas>, BidirectionalFMIndex<TSpec, TConfig> > & index,
                        StringSet<TText, TIrgendwas> const & text,
                        FibreSALF const & fibre,
                        TLambda const & progressCallback)
{
	// TODO:cpockrandt: progressCallback funzt so nicht, aber erstmal egal!
    //return indexCreate(index.fwd, text, fibre, progressCallback);
    return indexCreate(index.fwd, text, fibre, progressCallback) &&
    		indexCreate(index.rev, index.revText, fibre, progressCallback);
}

template <typename TText, typename TSpec, typename TLambda>
inline bool indexCreate(Index<TText, IndexSa<TSpec> > & index,
                        TText const & text,
                        FibreSA const &,
                        TLambda const & progressCallback)
{
    typedef Index<TText, IndexSa<TSpec> >      TIndex;
    typedef typename Fibre<TIndex, FibreSA>::Type       TSA;
    typedef typename DefaultIndexCreator<TIndex, FibreSA>::Type  TAlgo;

    indexText(index) = text;

    if (empty(text))
        return false;

    TSA & sa = getFibre(index, FibreSA());

    // Create the full SA.
    resize(sa, lengthSum(text), Exact());
    createSuffixArray(sa, text, TAlgo(), progressCallback);

    return true;
}

// ----------------------------------------------------------------------------
// Generic Sequence loading
// ----------------------------------------------------------------------------

template <typename TString, typename TSpec, typename TFormat>
inline int
loadSequences(StringSet<TString, TSpec > & seqs,
              CharString const & path,
              TFormat const & /*tag*/)
{
    //TODO experiment with differen file types
    std::ifstream stream;
    stream.open(toCString(path));
    if (!stream.is_open())
        return -1;

    typedef RecordReader<std::ifstream, DoublePass<> > TReader;
    TReader reader(stream);

    StringSet<CharString, TSpec > ids;

    int res = read2(ids, seqs, reader, TFormat());
    if (res)
        std::cerr << "Error : " << res << "\n";

    stream.close();
    return res;
}

template <typename TString, typename TSpec, typename TFormat>
inline int
loadIds(StringSet<TString, TSpec > & ids,
              CharString const & path,
              TFormat const & /*tag*/)
{
    //TODO experiment with differen file types
    std::ifstream stream;
    stream.open(toCString(path));
    if (!stream.is_open())
        return -1;

    typedef RecordReader<std::ifstream, DoublePass<> > TReader;
    TReader reader(stream);

    StringSet<CharString, TSpec > seqs;

    int res = read2(ids, seqs, reader, TFormat());
    if (res)
        std::cerr << "Error : " << res << "\n";

    stream.close();
    return res;
}

template <typename TSeqs, typename TIds, typename TFormat>
inline int
loadSeqsAndIds(TIds             & ids,
               TSeqs            & seqs,
               CharString const & path,
               TFormat    const & /*tag*/)
{
    //TODO experiment with differen file types
    std::ifstream stream;
    stream.open(toCString(path));
    if (!stream.is_open())
        return -1;

    typedef RecordReader<std::ifstream, DoublePass<> > TReader;
    TReader reader(stream);

    int res = read2(ids, seqs, reader, TFormat());
    if (res)
        std::cerr << "Error : " << res << "\n";

    stream.close();
    return res;
}


template <typename TString, typename TSpec, typename TFormat>
inline int
loadIds2(StringSet<TString, TSpec > & ids,
         CharString const & path,
         TFormat const & /*tag*/)
{
    //TODO experiment with differen file types
    std::ifstream stream;
    stream.open(toCString(path));
    if (!stream.is_open())
        return -1;

    typedef RecordReader<std::ifstream, SinglePass<> > TReader;
    TReader reader(stream);

    int res = 0;
    CharString seq;
    for (unsigned i = 0; i < length(ids); ++i)
    {
//         res = readRecord(ids[i], seq, reader, Fasta());
        if (value(reader) != '>')
            return 9;

        res = goNext(reader);
        if (res)
            std::cerr << "Error : " << res << "\n";
        res = skipBlanks(reader);
        if (res)
            std::cerr << "Error : " << res << "\n";
        res = readLine(ids[i], reader);
        if (res)
            std::cerr << "Error : " << res << "\n";
        res = skipUntilLineBeginsWithChar(reader, '>');
        if ((res) && (res != EOF_BEFORE_SUCCESS))
            std::cerr << "Error : " << res << "\n";
    }
    stream.close();
    return 0;
}


// ----------------------------------------------------------------------------
// truncate sequences
// ----------------------------------------------------------------------------

// template <typename TString, typename TSpec>
// inline void
// _debug_shorten(StringSet<TString, TSpec > & seqs, unsigned const len)
// {
//     StringSet<TString, TSpec > copySeqs;
//     reserve(copySeqs.concat, length(seqs)*len, Exact());
// 
//     for (TString const & s : seqs)
//         if (length(s) >= len)
//             appendValue(copySeqs, prefix(s, len), Exact());
// 
//     clear(seqs);
//     reserve(seqs.concat, length(copySeqs)*len, Exact());
//     for (TString const & s : copySeqs)
//         appendValue(seqs, s);
// }


// ----------------------------------------------------------------------------
// print if certain verbosity is set
// ----------------------------------------------------------------------------


template <typename T>
inline void
myPrintImpl(SharedOptions const & /**/,
            T const & first)
{
    std::cout << first;
}

inline void
myPrintImpl(SharedOptions const & options,
            std::stringstream const & first)
{
    std::string str = first.str();
//     std::cerr << "terminal cols: " << options.terminalCols
//               << " str.size() " << str.size() << "\n";
    if (options.isTerm && (str.size() >= (options.terminalCols -12)))
        std::cout << str.substr(str.size()-options.terminalCols+12,
                                options.terminalCols);
    else
        std::cout << str;
}

template <typename T, typename ... Args>
inline void
myPrintImpl(SharedOptions const & options,
            T const & first,
            Args const & ... args)
{
    myPrintImpl(options, first);
    myPrintImpl(options, args...);
}

template <typename ... Args>
inline void
myPrintImplThread(SharedOptions const & options,
//                   T const & first,
                  Args const & ... args)
{
    #pragma omp critical(stdout)
    {
//                 std::cout << "\033[" << omp_get_thread_num() << "B";
//                 std::cout << "\033E";
        if (options.isTerm)
        {
            for (unsigned char i=0; i< omp_get_thread_num(); ++i)
                std::cout << std::endl;
            std::cout << "\033[K";
        }
        std::cout << "Thread " << std::setw(3) << omp_get_thread_num() << "| ";

        myPrintImpl(options, args...);
        std::cout << "\n" << std::flush;
        if (options.isTerm)
            std::cout << "\033[" << omp_get_thread_num()+1 << "A";
    }
}

template <typename... Args>
inline void
myPrint(SharedOptions const & options, const int verbose, Args const &... args)
{
    if (options.verbosity >= verbose)
    {
        #if defined(_OPENMP)
        if (omp_in_parallel())
            myPrintImplThread(options, args...);
        else
        #endif
            myPrintImpl(options, args...);

        std::cout << std::flush;
    }
}

template <typename T>
inline void
appendToStatusImpl(std::stringstream & status,
                   T const & first)
{
    status << first;
}

template <typename T, typename ... Args>
inline void
appendToStatusImpl(std::stringstream & status,
                   T const & first,
                   Args const & ... args)
{
    appendToStatusImpl(status, first);
    appendToStatusImpl(status, args...);
}

template <typename... Args>
inline void
appendToStatus(std::stringstream & status,
               LambdaOptions const & options,
               const int verbose,
               Args const & ... args)
{
    if (options.verbosity >= verbose)
        appendToStatusImpl(status, args...);
}

// ----------------------------------------------------------------------------
// remove tag type
// ----------------------------------------------------------------------------

// template <typename T>
// T unTag(Tag<T> const & /**/)
// {
//     return T();
// }

// ----------------------------------------------------------------------------
// get plus-minus-range with bounds-checking for unsigned types
// ----------------------------------------------------------------------------

// template <typename TNum, typename TNum2>
// inline TNum
// _protectUnderflow(const TNum n, const TNum2 s)
// {
//     const TNum r = n -s;
//     return std::min(r, n);
// }
// 
// template <typename TNum, typename TNum2>
// inline TNum
// _protectOverflow(const TNum n, const TNum2 s)
// {
//     const TNum r = n + s;
//     return std::max(r, n);
// }
// 
// template <typename TGaps>
// inline bool
// _startsWithGap(TGaps const & gaps)
// {
//     SEQAN_ASSERT(length(gaps._array) > 0);
//     return (gaps._array[0] != 0);
// }
// 
// template <typename TGaps>
// inline int
// _endsWithGap(TGaps const & gaps)
// {
//     SEQAN_ASSERT(length(gaps._array) > 0);
//     if ((length(gaps._array)-1) % 2 == 1)
//         return -1;
//     return ((gaps._array[length(gaps._array)-1] != 0) ? 1 : 0);
// }
// 
// template <typename TGaps, typename TSeq>
// inline void
// _prependNonGaps(TGaps & gaps, TSeq const & seq)
// {
//     if (_startsWithGap(gaps))
//     {
//         insertValue(gaps._array, 0, length(seq)); // new non-gap column
//         insertValue(gaps._array, 0, 0); // empty gaps column
//     }
//     else
//     {
//         gaps._array[1] += length(seq);
//     }
// 
//     insert(value(gaps._source), 0, seq);
//     setBeginPosition(gaps, 0);
// }
// 
// template <typename TGaps, typename TSeq>
// inline void
// _appendNonGaps(TGaps & gaps, TSeq const & seq)
// {
//     switch (_endsWithGap(gaps))
//     {
//         case -1:
//         case  1:
//             appendValue(gaps._array, length(seq)); // new non-gap column
//             break;
//         case 0:
//             gaps._array[1] += length(seq);
//     }
//     append(value(gaps._source), seq);
//     setEndPosition(gaps, length(value(gaps._source)));
// }

#endif // header guard
