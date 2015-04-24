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
// lambda.cpp: Main File for the main application
// ==========================================================================

// why is this neccessary?
#undef SEQAN_HAS_ZLIB

// #define SEQAN_DEBUG_INDEX

// #define PARALLEL_SORT 0
// 0 = off
// 1 = GCC
// 2 = omptl

// #if PARALLEL_SORT == 1
//     #include <parallel/algorithm>
//     #define SORT __gnu_parallel::sort
// #elif PARALLEL_SORT == 2
//     #include <omptl/omptl_algorithm>
//     #define SORT omptl::sort
// #else
//     #define SORT std::sort
// #endif

#if defined(__GNUG__) && defined(_OPENMP)
    #define _GLIBCXX_PARALLEL
#endif
#define _GLIBCXX_USE_C99 1

#include <seqan/basic.h>

#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>

#include <iostream>

#include "lambda_indexer.hpp"

using namespace seqan;

// ==========================================================================
// Functions
// ==========================================================================
// --------------------------------------------------------------------------
// Function main()
// --------------------------------------------------------------------------


template <BlastFormatProgram p>
int mainTyped(LambdaIndexerOptions const & options,
              BlastFormat<BlastFormatFile::INVALID_File, p,
                BlastFormatGeneration::INVALID_Generation> const & /*tag*/);

template <BlastFormatProgram p,
          typename TRedAlph>
inline int
mainAlphed(TRedAlph const & /**/,
           LambdaIndexerOptions const & options,
           BlastFormat<BlastFormatFile::INVALID_File, p,
                BlastFormatGeneration::INVALID_Generation> const & /*tag*/);

template <BlastFormatProgram p,
          typename TRedAlph,
          typename TIndexSpecSpec>
inline int
mainIndexTyped(TRedAlph const & /**/,
               LambdaIndexerOptions const & options,
               TIndexSpecSpec const & /**/,
               BlastFormat<BlastFormatFile::INVALID_File, p,
                 BlastFormatGeneration::INVALID_Generation> const & /*tag*/);

// Program entry point.
int main(int argc, char const ** argv)
{
    // Parse the command line.
    seqan::ArgumentParser parser;
    LambdaIndexerOptions options;
    seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

    // If there was an error parsing or built-in argument parser functionality
    // was triggered then we exit the program.  The return code is 1 if there
    // were errors and 0 if there were none.
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;

  // CONVERT Run-time options to compile-time Format-Type
    switch (options.blastProg)
    {
        case BlastFormatProgram::BLASTN :
            {
                typedef BlastFormat<BlastFormatFile::INVALID_File,
                                    BlastFormatProgram::BLASTN,
                                    BlastFormatGeneration::INVALID_Generation> format;
                return mainTyped(options, format());
            }break;
        case BlastFormatProgram::BLASTP :
            {
                typedef BlastFormat<BlastFormatFile::INVALID_File,
                                    BlastFormatProgram::BLASTP,
                                    BlastFormatGeneration::INVALID_Generation> format;
                return mainTyped(options, format());
            }
            break;
        case BlastFormatProgram::BLASTX :
            {
                typedef BlastFormat<BlastFormatFile::INVALID_File,
                                    BlastFormatProgram::BLASTX,
                                    BlastFormatGeneration::INVALID_Generation> format;
                return mainTyped(options, format());
            }
            break;
        case BlastFormatProgram::TBLASTN :
            {
                typedef BlastFormat<BlastFormatFile::INVALID_File,
                                    BlastFormatProgram::TBLASTN,
                                    BlastFormatGeneration::INVALID_Generation> format;
                return mainTyped(options, format());
            }
            break;
        case BlastFormatProgram::TBLASTX :
            {
                typedef BlastFormat<BlastFormatFile::INVALID_File,
                                    BlastFormatProgram::TBLASTX,
                                    BlastFormatGeneration::INVALID_Generation> format;
                return mainTyped(options, format());
            }
            break;
        default:
            return -1;
    }
    return -1;
}


template <BlastFormatProgram p>
int mainTyped(LambdaIndexerOptions const & options,
              BlastFormat<BlastFormatFile::INVALID_File, p,
                BlastFormatGeneration::INVALID_Generation> const & /*tag*/)
{
    using TFormat   = BlastFormat<BlastFormatFile::INVALID_File,
                                  p,
                                  BlastFormatGeneration::INVALID_Generation>;
    myPrint(options, 1, "LAMBDA - Indexer"
                      "\n======================================================"
                      "\nVersion ", LAMBDA_VERSION, "\n\n");

    switch (options.alphReduction)
    {
        case 0:
        {
            typedef AminoAcid TAA;
            return mainAlphed(TAA(), options, TFormat());
        } break;
//         case 1:
//         {
//             typedef AminoAcid10 TAA;
//             return mainAlphed(TAA(), options, TFormat());
//         } break;
        case 2:
        {
            typedef ReducedAminoAcid<Murphy10> TAA;
            return mainAlphed(TAA(), options, TFormat());
        } break;
//         case 8:
//         {
//             typedef ReducedAminoAcid<ClusterReduction<8>> TAA;
//             return mainAlphed(TAA(), options, TFormat());
//         } break;
// 
//         case 10:
//         {
//             typedef ReducedAminoAcid<ClusterReduction<10>> TAA;
//             return mainAlphed(TAA(), options, TFormat());
//         } break;
// 
//         case 12:
//         {
//             typedef ReducedAminoAcid<ClusterReduction<12>> TAA;
//             return mainAlphed(TAA(), options, TFormat());
//         } break;
    }
    return -1;
}

template <BlastFormatProgram p,
          typename TRedAlph>
inline int
mainAlphed(TRedAlph const & /**/,
           LambdaIndexerOptions const & options,
           BlastFormat<BlastFormatFile::INVALID_File, p,
                BlastFormatGeneration::INVALID_Generation> const & /*tag*/)
{
    using TFormat   = BlastFormat<BlastFormatFile::INVALID_File,
                                  p,
                                  BlastFormatGeneration::INVALID_Generation>;

    if (options.algo == "mergesort")
        return mainIndexTyped(TRedAlph(), options, SaAdvancedSort<MergeSortTag>(), TFormat());
    else if (options.algo == "quicksort")
        return mainIndexTyped(TRedAlph(), options, SaAdvancedSort<QuickSortTag>(), TFormat());
    else if (options.algo == "quicksortbuckets")
        return mainIndexTyped(TRedAlph(), options, SaAdvancedSort<QuickSortBucketTag>(), TFormat());
    else
        return mainIndexTyped(TRedAlph(), options, Nothing(), TFormat());
}

template <BlastFormatProgram p,
          typename TRedAlph,
          typename TIndexSpecSpec>
inline int
mainIndexTyped(TRedAlph const & /**/,
               LambdaIndexerOptions const & options,
               TIndexSpecSpec const & /**/,
               BlastFormat<BlastFormatFile::INVALID_File, p,
                 BlastFormatGeneration::INVALID_Generation> const & /*tag*/)
{
    using TFormat   = BlastFormat<BlastFormatFile::INVALID_File,
                                  p,
                                  BlastFormatGeneration::INVALID_Generation>;

    using TOrigSet  = TCDStringSet<OrigSubjAlph<p>>;
    using TTransSet = TCDStringSet<TransAlph<p>>;

    TTransSet translatedSeqs;

    {
        TOrigSet originalSeqs;
        int ret = 0;

        // ids get saved to disk again immediately and are not kept in memory
        ret = loadSubjSeqsAndIds(originalSeqs, options);
        if (ret)
            return ret;

        // preserve lengths of untranslated sequences
        _saveOriginalSeqLengths(originalSeqs.limits,
                                options,
                                SIsTranslated<TFormat>());

        // convert the seg file to seqan binary format
        ret = convertMaskingFile(length(originalSeqs), options);
        if (ret)
            return ret;

        // translate or swap depending on program
        translateOrSwap(translatedSeqs, originalSeqs, options);
    }

    // dump translated and unreduced sequences
    dumpTranslatedSeqs(translatedSeqs, options);

    // see if final sequence set actually fits into index 
    if (!checkIndexSize(translatedSeqs))
        return -1;

    if (options.dbIndexType == 2)
    {
        using TIndexSpec = TBidirectionalFMIndex<TIndexSpecSpec>;
        generateIndexAndDump<TIndexSpec,TIndexSpecSpec>(translatedSeqs,
                                                        options,
                                                        TRedAlph(),
                                                        TFormat());
    } else if (options.dbIndexType == 1)
	{
		using TIndexSpec = TFMIndex<TIndexSpecSpec>;
		generateIndexAndDump<TIndexSpec,TIndexSpecSpec>(translatedSeqs,
														options,
														TRedAlph(),
														TFormat());
	} else
    {
        using TIndexSpec = IndexSa<TIndexSpecSpec>;
        generateIndexAndDump<TIndexSpec,TIndexSpecSpec>(translatedSeqs,
                                                        options,
                                                        TRedAlph(),
                                                        TFormat());
    }

    return 0;
}

