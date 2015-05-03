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
// lambda.cpp: Main File for Lambda
// ==========================================================================

#define FASTBUILD

#define LOOKUP_TABLE_SIZE	24*24*24*4

// #define SEQAN_DEBUG_INDEX

#define _GLIBCXX_USE_C99 1
//#define LAMBDA_BITCOPMRESSED_STRINGS 1

#include <vector>
#include <iostream>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>
#include <seqan/reduced_aminoacid.h>
#include <seqan/misc/misc_terminal.h>


#include <seqan/index.h>
//#include <seqan/index/find_index_bidirectional.h>

#include "options.hpp"
#include "match.hpp"
#include "lambda.hpp"
#include "misc.hpp"

using namespace seqan;


inline BlastFormatFile
_fileType(LambdaOptions const & options)
{
    if (endsWith(options.output, ".m0"))
        return BlastFormatFile::PAIRWISE;
    else if (endsWith(options.output, ".m8"))
        return BlastFormatFile::TABULAR;
    else if (endsWith(options.output, ".m9"))
        return BlastFormatFile::TABULAR_WITH_HEADER;
    else
        return BlastFormatFile::INVALID_File;
}

// forwards

inline int
argConv0(LambdaOptions const & options);
//-
template <BlastFormatFile m,
          BlastFormatGeneration g>
inline int
argConv1(LambdaOptions                               const & options,
         BlastFormat<m,BlastFormatProgram::BLASTN,g> const & /**/);
template <BlastFormatFile m,
          BlastFormatProgram p,
          BlastFormatGeneration g,
          MyEnableIf<p != BlastFormatProgram::BLASTN> = 0>
inline int
argConv1(LambdaOptions      const & options,
         BlastFormat<m,p,g> const & /**/);
//-
template <BlastFormatFile m,
          BlastFormatProgram p,
          BlastFormatGeneration g,
          typename TRedAlph>
inline int
argConv2(LambdaOptions      const & options,
         BlastFormat<m,p,g> const & /**/,
         TRedAlph           const & /**/);
//-
template <BlastFormatFile m,
          BlastFormatProgram p,
          BlastFormatGeneration g,
          typename TRedAlph,
          typename TScoreScheme>
inline int
argConv3(LambdaOptions      const & options,
         BlastFormat<m,p,g> const & /**/,
         TRedAlph           const & /**/,
         TScoreScheme       const & /**/);
//-
template <BlastFormatFile m,
          BlastFormatProgram p,
          BlastFormatGeneration g,
          typename TRedAlph,
          typename TScoreScheme,
          typename TScoreExtension>
inline int
preMain(LambdaOptions      const & options,
         BlastFormat<m,p,g> const & /**/,
         TRedAlph           const & /**/,
         TScoreScheme       const & /**/,
         TScoreExtension    const & /**/);
//-
template <typename TIndexSpec,
          typename TRedAlph,
          typename TScoreScheme,
          typename TScoreExtension,
          BlastFormatFile m,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline int
realMain(LambdaOptions      const & options,
         BlastFormat<m,p,g> const & /**/,
         TRedAlph           const & /**/,
         TScoreScheme       const & /**/,
         TScoreExtension    const & /**/);

// --------------------------------------------------------------------------
// Function main()
// --------------------------------------------------------------------------

template <typename TAlphabet>
void showAllLettersOfMyAlphabet(TAlphabet const &)
{
    typedef typename Size<TAlphabet>::Type TSize;
    TSize alphSize = ValueSize<TAlphabet>::VALUE;
    for (TSize i = 0; i < alphSize; ++i)
        std::cout << i << ',' << TAlphabet(i) << "  ";
    std::cout << std::endl;
}


// Program entry point.

int main(int argc, char const ** argv)
{
	//measureCPUTime();
	//return 0;
	/*for (int n = 0; n < 20; ++n)
	{
		std::cout << "n = " << n << " ";
		for (int p = 0; p <= n; ++p)
		{
			std::cout << ".";
			for (int i = 0; i < 100; ++i)
			{
				if (testExactSearch(n, p, false))
					return 1;
				if (testExactSearch2(n, p, 1, false))
					return 1;
				if (testExactSearch2(n, p, 2, false))
					return 1;
				if (testExactSearch2(n, p, 5, false))
					return 1;
			}
		}
		std::cout << std::endl;
	}
	std::cout << "n = 10000, p = 30" << std::endl;
	for (int i = 0; i < 100; ++i)
	{
		std::cout << i << std::endl;
		if (testExactSearch(10000, 50, false))
			break;
		if (testExactSearch2(10000, 50, 100, false))
			break;
	}

	return 0;*/

    // Parse the command line.
    seqan::ArgumentParser parser;
    LambdaOptions options;
    seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

    // If there was an error parsing or built-in argument parser functionality
    // was triggered then we exit the program.  The return code is 1 if there
    // were errors and 0 if there were none.
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;

//     std::cout <<   "Match  sizeof : " << sizeof(Match)
//               << "\n       alignof: " << alignof(Match)
//               << "\n    is_trivial: " << std::is_trivial<Match>::value
// //               << "\ntrivially_copy: " << std::is_trivially_copyable<Match>::value
//               << "\n";

    return argConv0(options);
}


// CONVERT Run-time options to compile-time Format-Type
inline int
argConv0(LambdaOptions const & options)
{
    switch (options.blastProg)
    {
         case BlastFormatProgram::BLASTN :
         {
 //             template <BlastFormatFile m>
 //             using TBF = BlastFormat<m,
 //                                 BlastFormatProgram::BLASTN,
 //                                 BlastFormatGeneration::BLAST>;
             switch (_fileType(options))
             {
                 case BlastFormatFile::PAIRWISE:
                 {
                     typedef BlastFormat<BlastFormatFile::PAIRWISE,
                                         BlastFormatProgram::BLASTN,
                                         BlastFormatGeneration::BLAST> TFormat;
                     return argConv1(options, TFormat());
                 } break;
                 case BlastFormatFile::TABULAR:
                 {
                     typedef BlastFormat<BlastFormatFile::TABULAR,
                                         BlastFormatProgram::BLASTN,
                                         BlastFormatGeneration::BLAST> TFormat;
                     return argConv1(options, TFormat());
                 } break;
                 case BlastFormatFile::TABULAR_WITH_HEADER:
                 {
                     typedef BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                                         BlastFormatProgram::BLASTN,
                                         BlastFormatGeneration::BLAST> TFormat;
                     return argConv1(options, TFormat());
                 } break;
                 default:
                     break;
             }
         } break;
#ifndef FASTBUILD
        case BlastFormatProgram::BLASTP :
        {
//             template <BlastFormatFile m>
//             using TBF = BlastFormat<m,
//                                 BlastFormatProgram::BLASTP,
//                                 BlastFormatGeneration::BLAST>;
            switch (_fileType(options))
            {
                case BlastFormatFile::PAIRWISE:
                {
                    typedef BlastFormat<BlastFormatFile::PAIRWISE,
                                        BlastFormatProgram::BLASTP,
                                        BlastFormatGeneration::BLAST> TFormat;
                    return argConv1(options, TFormat());
                } break;
                case BlastFormatFile::TABULAR:
                {
                    typedef BlastFormat<BlastFormatFile::TABULAR,
                                        BlastFormatProgram::BLASTP,
                                        BlastFormatGeneration::BLAST> TFormat;
                    return argConv1(options, TFormat());
                } break;
                case BlastFormatFile::TABULAR_WITH_HEADER:
                {
                    typedef BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                                        BlastFormatProgram::BLASTP,
                                        BlastFormatGeneration::BLAST> TFormat;
                    return argConv1(options, TFormat());
                } break;
                default:
                    break;
            }
        } break;
#endif
        case BlastFormatProgram::BLASTX :
        {
//             template <BlastFormatFile m>
//             using TBF = BlastFormat<m,
//                                 BlastFormatProgram::BLASTX,
//                                 BlastFormatGeneration::BLAST>;
            switch (_fileType(options))
            {
// #ifndef FASTBUILD
                case BlastFormatFile::PAIRWISE:
                {
                    typedef BlastFormat<BlastFormatFile::PAIRWISE,
                                        BlastFormatProgram::BLASTX,
                                        BlastFormatGeneration::BLAST> TFormat;
                    return argConv1(options, TFormat());
                } break;
// #endif
                case BlastFormatFile::TABULAR:
                {
                    typedef BlastFormat<BlastFormatFile::TABULAR,
                                        BlastFormatProgram::BLASTX,
                                        BlastFormatGeneration::BLAST> TFormat;
                    return argConv1(options, TFormat());
                } break;
// #ifndef FASTBUILD
                case BlastFormatFile::TABULAR_WITH_HEADER:
                {
                    typedef BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                                        BlastFormatProgram::BLASTX,
                                        BlastFormatGeneration::BLAST> TFormat;
                    return argConv1(options, TFormat());
                } break;
// #endif
                default:
                    break;
            }
        } break;
#ifndef FASTBUILD
        case BlastFormatProgram::TBLASTN :
        {
//             template <BlastFormatFile m>
//             using TBF = BlastFormat<m,
//                                 BlastFormatProgram::TBLASTN,
//                                 BlastFormatGeneration::BLAST>;
            switch (_fileType(options))
            {
                case BlastFormatFile::PAIRWISE:
                {
                    typedef BlastFormat<BlastFormatFile::PAIRWISE,
                                        BlastFormatProgram::TBLASTN,
                                        BlastFormatGeneration::BLAST> TFormat;
                    return argConv1(options, TFormat());
                } break;
                case BlastFormatFile::TABULAR:
                {
                    typedef BlastFormat<BlastFormatFile::TABULAR,
                                        BlastFormatProgram::TBLASTN,
                                        BlastFormatGeneration::BLAST> TFormat;
                    return argConv1(options, TFormat());
                } break;
                case BlastFormatFile::TABULAR_WITH_HEADER:
                {
                    typedef BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                                        BlastFormatProgram::TBLASTN,
                                        BlastFormatGeneration::BLAST> TFormat;
                    return argConv1(options, TFormat());
                } break;
                default:
                    break;
            }
        } break;

        case BlastFormatProgram::TBLASTX :
        {
//             template <BlastFormatFile m>
//             using TBF = BlastFormat<m,
//                                 BlastFormatProgram::TBLASTX,
//                                 BlastFormatGeneration::BLAST>;
            switch (_fileType(options))
            {
                case BlastFormatFile::PAIRWISE:
                {
                    typedef BlastFormat<BlastFormatFile::PAIRWISE,
                                        BlastFormatProgram::TBLASTX,
                                        BlastFormatGeneration::BLAST> TFormat;
                    return argConv1(options, TFormat());
                } break;
                case BlastFormatFile::TABULAR:
                {
                    typedef BlastFormat<BlastFormatFile::TABULAR,
                                        BlastFormatProgram::TBLASTX,
                                        BlastFormatGeneration::BLAST> TFormat;
                    return argConv1(options, TFormat());
                } break;
                case BlastFormatFile::TABULAR_WITH_HEADER:
                {
                    typedef BlastFormat<BlastFormatFile::TABULAR_WITH_HEADER,
                                        BlastFormatProgram::TBLASTX,
                                        BlastFormatGeneration::BLAST> TFormat;
                    return argConv1(options, TFormat());
                } break;
                default:
                    break;
            }
        } break;
#endif
        default:
            break;
    }
    return -1;
}

/// Alphabet reduction

template <BlastFormatFile m,
          BlastFormatGeneration g>
inline int
argConv1(LambdaOptions                               const & options,
         BlastFormat<m,BlastFormatProgram::BLASTN,g> const & /**/)
{
    using TFormat = BlastFormat<m,BlastFormatProgram::BLASTN,g>;
    return argConv2(options, TFormat(), Dna5());
}

template <BlastFormatFile m,
          BlastFormatProgram p,
          BlastFormatGeneration g,
          MyEnableIf<p != BlastFormatProgram::BLASTN> >
inline int
argConv1(LambdaOptions      const & options,
         BlastFormat<m,p,g> const & /**/)
{
    using TFormat = BlastFormat<m,p,g>;
    switch (options.alphReduction)
    {

        case 0:
            return argConv2(options, TFormat(), AminoAcid());
        case 2:
            return argConv2(options, TFormat(), ReducedAminoAcid<Murphy10>());
#if 0
        case 10:
            return argConv2(options, TFormat(), ReducedAminoAcid<ClusterReduction<10>>());
        case 1:
            return argConv2(options, TFormat(), AminoAcid10());
        case 8:
            return argConv2(options, TFormat(), ReducedAminoAcid<ClusterReduction<8>>());
        case 12:
            return argConv2(options, TFormat(), ReducedAminoAcid<ClusterReduction<12>>());
#endif
        default:
            return -1;
    }
    return -1;
}

/// scoring scheme

template <BlastFormatFile m,
          BlastFormatProgram p,
          BlastFormatGeneration g,
          typename TRedAlph>
inline int
argConv2(LambdaOptions      const & options,
         BlastFormat<m,p,g> const & /**/,
         TRedAlph           const & /**/)
{
    using TFormat = BlastFormat<m,p,g>;
    switch (options.scoringMethod)
    {

        case 0:
            return argConv3(options, TFormat(), TRedAlph(), Score<int, Simple>());
#ifndef FASTBUILD
        case 45:
            return argConv3(options, TFormat(), TRedAlph(), Blosum45());
        case 80:
            return argConv3(options, TFormat(), TRedAlph(), Blosum80());
#endif
        case 62:
            return argConv3(options, TFormat(), TRedAlph(), Blosum62());
        default:
            std::cerr << "Unsupported Scoring Scheme selected.\n";
            break;
    }

    return -1;
}

template <BlastFormatFile m,
          BlastFormatProgram p,
          BlastFormatGeneration g,
          typename TRedAlph,
          typename TScoreScheme>
inline int
argConv3(LambdaOptions      const & options,
         BlastFormat<m,p,g> const & /**/,
         TRedAlph           const & /**/,
         TScoreScheme       const & /**/)
{
    using TFormat = BlastFormat<m,p,g>;
#ifndef FASTBUILD
    if (options.gapOpen == 0)
        return preMain(options,
                        TFormat(),
                        TRedAlph(),
                        TScoreScheme(),
                        LinearGaps());
    else
#endif
        return preMain(options,
                        TFormat(),
                        TRedAlph(),
                        TScoreScheme(),
                        AffineGaps());

}

template <BlastFormatFile m,
          BlastFormatProgram p,
          BlastFormatGeneration g,
          typename TRedAlph,
          typename TScoreScheme,
          typename TScoreExtension>
inline int
preMain(LambdaOptions      const & options,
         BlastFormat<m,p,g> const & /**/,
         TRedAlph           const & /**/,
         TScoreScheme       const & /**/,
         TScoreExtension    const & /**/)
{
    using TFormat = BlastFormat<m,p,g>;
    int indexType = options.dbIndexType;
//     if (indexType == -1) // autodetect
//     {
//         //TODO FIX THIS WITH NEW EXTENSIONS
//         CharString file = options.dbFile;
//         append(file, ".sa");
//         struct stat buffer;
//         if (stat(toCString(file), &buffer) == 0)
//         {
//             indexType = 0;
//         } else
//         {
//             file = options.dbFile;
//             append(file, ".sa.val"); // FM Index
//             struct stat buffer;
//             if (stat(toCString(file), &buffer) == 0)
//             {
//                 indexType = 1;
//             } else
//             {
//                 std::cerr << "No Index file could be found, please make sure paths "
//                         << "are correct and the files are readable.\n" << std::flush;
// 
//                 return -1;
//             }
//         }
//     }

    /*if (indexType == 0)
        return realMain<IndexSa<>>(options,
                                   TFormat(),
                                   TRedAlph(),
                                   TScoreScheme(),
                                   TScoreExtension());
    else if (indexType == 1)
        return realMain<TFMIndex<>>(options,
                                   TFormat(),
                                   TRedAlph(),
                                   TScoreScheme(),
                                   TScoreExtension());
    else*/
    if (indexType == 2)
        return realMain<TBidirectionalFMIndex<>>(options,
                                   TFormat(),
                                   AminoAcid(),/*TRedAlph(),*/
                                   TScoreScheme(),
                                   TScoreExtension());
    std::cout << "aahhhhhh! " << indexType << std::endl;
    return 3;
}

/// REAL MAIN
#ifdef _OPENMP
#define TID omp_get_thread_num()
#else
#define TID 0
#endif
template <typename TIndexSpec,
          typename TRedAlph,
          typename TScoreScheme,
          typename TScoreExtension,
          BlastFormatFile m,
          BlastFormatProgram p,
          BlastFormatGeneration g>
inline int
realMain(LambdaOptions      const & options,
         BlastFormat<m,p,g> const & /**/,
         TRedAlph           const & /**/,
         TScoreScheme       const & /**/,
         TScoreExtension    const & /**/)
{
    using TGlobalHolder = GlobalDataHolder<TRedAlph,
                                           TScoreScheme,
                                           TIndexSpec,
                                           m, p, g>;
    using TLocalHolder = LocalDataHolder<Match, TGlobalHolder, TScoreExtension>;






    myPrint(options, 1, "LAMBDA - the Local Aligner for Massive Biological DatA"
                      "\n======================================================"
                      "\nVersion ", LAMBDA_VERSION, "\n\n");

    if (options.verbosity >= 2)
        printOptions<TLocalHolder>(options);

    TGlobalHolder globalHolder;

    /*typedef StringSet<String<SimpleType<unsigned char, Dna5_>, Alloc<> >, Owner<ConcatDirect<void> > > TMyText;
    typedef BidirectionalFMIndex<void, FMIndexConfig<void, FMBidirectional> > TMyIndex;
    //Index<TMyText, TMyIndex> testIndex("aa");
    //indexSA(testIndex);
    return 0;*/

    int ret = prepareScoring(globalHolder, options);
    if (ret)
        return ret;

    ret = loadSubjects(globalHolder, options);
    if (ret)
        return ret;

    ret = loadDbIndexFromDisk(globalHolder, options);
    if (ret)
        return ret;

    ret = loadSegintervals(globalHolder, options);
    if (ret)
        return ret;

    ret = loadQuery(globalHolder, options);
    if (ret)
        return ret;

    std::ofstream stream;
    stream.open(toCString(options.output));
    if (!stream.is_open())
    {
        std::cerr << "Error opening output file for writing." << std::endl;
        return -2;
    }

    ret = writeTop(stream,
                   globalHolder.dbSpecs,
                   typename TGlobalHolder::TFormat());
    if (ret)
        return ret;



    unsigned long start3 = std::clock();


    /*unsigned long countAminoAcids[24];
    for (int z = 0; z < 24; ++z)
    	countAminoAcids[z] = 0;*/

    //std::vector<std::vector<unsigned long> > preComputedIndexPos(LOOKUP_TABLE_SIZE);
    unsigned long preComputedIndexPos[LOOKUP_TABLE_SIZE] = { 0 };

    typedef typename Size<TRedAlph>::Type TSize;
	TSize alphSize = ValueSize<TRedAlph>::VALUE;

    typedef typename Iterator<decltype(globalHolder.dbIndex), TopDown<> >::Type TIndexIt;
	TIndexIt indexIt(globalHolder.dbIndex);

	//unsigned int pos = 0;
    //for (TSize a = 0; a < alphSize; ++a)
    //{
    	//auto vDescA1 = value(indexIt.fwdIter);
    	//auto vDescA2 = value(indexIt.bwdIter);
    	//goDown(indexIt.fwdIter, TRedAlph(a));

	for (TSize b = 0; b < alphSize; ++b)
	{
		if (goDown(indexIt.fwdIter, TRedAlph(b)))
		{
			auto vDescB1 = value(indexIt.fwdIter);
			auto vDescB2 = value(indexIt.bwdIter);
			for (TSize c = 0; c < alphSize; ++c)
			{
				if (goDown(indexIt.fwdIter, TRedAlph(c)))
				{
					auto vDescC1 = value(indexIt.fwdIter);
					auto vDescC2 = value(indexIt.bwdIter);
					for (TSize d = 0; d < alphSize; ++d)
					{
						if (goDown(indexIt.fwdIter, TRedAlph(d)))
						{
							unsigned int pos = 4*b + 4*24*c + 4*24*24*d;
							preComputedIndexPos[pos] = indexIt.fwdIter.vDesc.range.i1;
							preComputedIndexPos[pos+1] = indexIt.fwdIter.vDesc.range.i2;
							preComputedIndexPos[pos+2] = indexIt.bwdIter.vDesc.range.i1;
							preComputedIndexPos[pos+3] = indexIt.bwdIter.vDesc.range.i2;
							value(indexIt.fwdIter) = vDescC1;
							value(indexIt.bwdIter) = vDescC2;
						}
						//++pos;
					}
					value(indexIt.fwdIter) = vDescB1;
					value(indexIt.bwdIter) = vDescB2;
				}
				else
				{
					//pos += alphSize;
				}
			}
		}
		else
		{
			//pos += alphSize*alpSize;
		}
		goRoot(indexIt);
	}
        	//value(indexIt.fwdIter) = vDescA1;
        	//value(indexIt.fwdIter) = vDescA2;
        //}
        //goRoot(indexIt);
    //}
    std::cout << "Clock-Ticks: " << (std::clock() - start3) << " ... " << ((std::clock() - start3)/CLOCKS_PER_SEC) << std::endl;









    if (options.doubleIndexing)
    {
        myPrint(options, 1,
                "Searching ",
                options.queryPart,
                " blocks of query with ",
                options.threads,
                " threads...\n");
        if ((options.isTerm) && (options.verbosity >= 1))
        {
            for (unsigned char i=0; i< options.threads+3; ++i)
                std::cout << std::endl;
            std::cout << "\033[" << options.threads+2 << "A";
        }
    } else
    {
        myPrint(options, 1, "Searching and extending hits on-line...progress:\n"
                "0%  10%  20%  30%  40%  50%  60%  70%  80%  90%  100%\n|");
    }
    double start = sysTime();

    // at least a block for each thread on double-indexing,
    // otherwise a block for each original query (contains 6 queries if
    // translation is used)
    uint64_t nBlocks = (options.doubleIndexing
                        ? options.queryPart
                        : length(globalHolder.qryIds));

    uint64_t lastPercent = 0;

    unsigned long start2 = std::clock();


    SEQAN_OMP_PRAGMA(parallel)
    {
        TLocalHolder localHolder(options, globalHolder);

        SEQAN_OMP_PRAGMA(for schedule(dynamic))
        for (uint64_t t = 0; t < nBlocks; ++t)
        {
            int res = 0;

            localHolder.init(t);

            // seed
            res = generateSeeds(localHolder/*, countAminoAcids*/);
            if (res)
                continue;

            if (options.doubleIndexing)
            {
                res = generateTrieOverSeeds(localHolder);
                if (res)
                    continue;
            }

            // search

            //double mystart, c1, c3;

            //mystart = std::clock();
            search(localHolder, preComputedIndexPos);
            //c1 = std::clock() - mystart;
            //std::cout << "\t\t\t\ta" << std::endl;

            // sort
            sortMatches(localHolder);
            //std::cout << "\t\t\t\tb" << std::endl;

            // extend
            //mystart = std::clock();
            res = iterateMatches(stream, localHolder);
            //c3 = std::clock() - mystart;
            //std::cout << "\t\t\t\tc" << std::endl;

            //std::cout << "" << c1 << "\t" << c3 << std::endl;
            if (res)
                continue;


            if ((!options.doubleIndexing) && (TID == 0) &&
                (options.verbosity >= 1))
            {
                unsigned curPercent = ((t * 50) / nBlocks) * 2; // round to even
                printProgressBar(lastPercent, curPercent);
            }

        } // implicit thread sync here

        if ((!options.doubleIndexing) && (TID == 0) && (options.verbosity >= 1))
            printProgressBar(lastPercent, 100);

        SEQAN_OMP_PRAGMA(critical(statsAdd))
        {
            globalHolder.stats += localHolder.stats;
        }
    }

    std::cout << "Clock-Ticks: " << (std::clock() - start2) << " ... " << ((std::clock() - start2)/CLOCKS_PER_SEC) << std::endl;

    /*std::cout << "-------------------" << std::endl;
    unsigned long totalAA = 0;
    for (int z = 0; z < 24; ++z)
    {
    	std::cout << AminoAcid(z) << ": " << countAminoAcids[z] << std::endl;
    	totalAA += countAminoAcids[z];
    }
    std::cout << "-------------------" << std::endl;
    std::cout << "Total: " << totalAA << std::endl;
    std::cout << "-------------------" << std::endl;*/

    ret = writeBottom(stream,
                      globalHolder.dbSpecs,
                      globalHolder.blastScoringAdapter,
                      typename TGlobalHolder::TFormat());

    stream.close();

    if (ret)
        return ret;

    if (!options.doubleIndexing)
    {
        myPrint(options, 2, "Runtime: ", sysTime() - start, "s.\n\n");
    }

    printStats(globalHolder.stats, options);

    return 0;
}
