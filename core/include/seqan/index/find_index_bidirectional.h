/*
 * find_index_todo.h
 *
 *  Created on: 25.01.2015
 *      Author: chris
 */

#include <seqan/seq_io.h>

#ifndef FIND_INDEX_TODO_H_
#define FIND_INDEX_TODO_H_

namespace seqan {

const std::string currentDateTime() {
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    strftime(buf, sizeof(buf), "%X", &tstruct);
    return buf;
}

template < typename TText, class TOccSpec, class TIndexSpec, typename TSpec >
inline void
_approxExtend(
	Iter<Index<TText, BidirectionalFMIndex<TOccSpec, TIndexSpec> >, VSTree<TopDown<ParentLinks<TSpec> > > > &it,
	int errors, TText const &str, std::list<Pair<int> > &ranges, bool leftExtension, unsigned int strPos = 0)
{
	typedef typename Value<TText>::Type						TChar;
	typedef Index<TText, FMIndex<TOccSpec, TIndexSpec> >	TIndex;
	typedef typename Value<TIndex>::Type                	TAlphabet;

	if (errors < 0 || strPos > length(str))
		return;
	if (strPos == length(str))
	{
		int i1 = value(it.fwdIter).range.i1;
		int i2 = value(it.fwdIter).range.i2;

		// merge adjacent intervals. Only need to check last element in list because search is performed in lexicographic order.
		if (ranges.back().i2 == i1)
			ranges.back().i2 = i2;
		else
			ranges.push_back(*(new Pair<int>(i1, i2)));

		return;
	}

	int alpSize = ValueSize<TAlphabet>::VALUE;
	for (TAlphabet i = 0; i < alpSize; ++i)
	{
		TChar c = TAlphabet(i);
		if (leftExtension && !leftExtend(it, c))
			continue;
		if (!leftExtension && !rightExtend(it, c))
			continue;

		if (c == str[strPos])
			_approxExtend(it, errors, str, ranges, leftExtension, strPos + 1);
		else
			_approxExtend(it, errors-1, str, ranges, leftExtension, strPos + 1);

		goUp(it.fwdIter);
	}
}

template < typename TText, class TOccSpec, class TIndexSpec, typename TSpec >
inline void
rightApproxExtend(
	Iter<Index<TText, BidirectionalFMIndex<TOccSpec, TIndexSpec> >, VSTree<TopDown<ParentLinks<TSpec> > > > &it,
	int errors, TText const &str, std::list<Pair<int> > &ranges)
{
	_approxExtend(it, errors, str, ranges, false, 0);
}

template < typename TText, class TOccSpec, class TIndexSpec, typename TSpec >
inline void
leftApproxExtend(
	Iter<Index<TText, BidirectionalFMIndex<TOccSpec, TIndexSpec> >, VSTree<TopDown<ParentLinks<TSpec> > > > &it,
	int errors, TText const &str, std::list<Pair<int> > &ranges)
{
	_approxExtend(it, errors, str, ranges, true, 0);
}

inline void genRandomStr(String<char> &s, const unsigned int len, const unsigned int alpSize)
{
	resize(s, len);
	for (unsigned int i = 0; i < len; ++i)
         s[i] = 'a' + rand()%alpSize;
}

inline bool
testExactSearch(unsigned int textLength, unsigned int patternLength, unsigned int alpSize, bool const debug = 0)
{
	typedef String<char>											TText;

	typedef Index<TText, FMIndex<> >								TFMIndex;
	typedef Iterator<TFMIndex, TopDown<> >::Type					TFMIter;

	typedef Index<TText, BidirectionalFMIndex<> >					TBiFMIndex;
	typedef Iterator<TBiFMIndex, TopDown<ParentLinks<> > >::Type	TBiFMIter;

	const char *fileName = "/home/chris/indexTest/myindex";

	unsigned int i;
	String<char> text;
	String<char> pattern;
	String<char> revText;

	srand(time(NULL));

	genRandomStr(text, textLength, alpSize);
	genRandomStr(pattern, patternLength, alpSize);

	ModifiedString< String<char>, ModReverse > revPattern(pattern);

	resize(revText, textLength);
	for (i = 0; i < textLength; ++i)
	{
		revText[i] = text[textLength - i - 1];
	}

	TFMIndex fmIndex1(text);
	TFMIndex fmIndex2(revText);
	TFMIter fm1(fmIndex1);
	TFMIter fm2(fmIndex2);

	TBiFMIndex bifmIndexOpen1(text);
	TBiFMIndex bifmIndexOpen2(text);

	indexRequire(bifmIndexOpen1.fwd, FibreSA());
	indexRequire(bifmIndexOpen1.rev, FibreSA());

	indexRequire(bifmIndexOpen2.fwd, FibreSA());
	indexRequire(bifmIndexOpen2.rev, FibreSA());

	system("exec rm -r /home/chris/indexTest/*");

	save(bifmIndexOpen1, fileName);
	save(bifmIndexOpen2, fileName);

	TBiFMIndex bifmIndex1(text);
	TBiFMIndex bifmIndex2(text);

	open(bifmIndex1, fileName);
	open(bifmIndex2, fileName);

	TBiFMIter bifm1(bifmIndex1);
	TBiFMIter bifm2(bifmIndex2);

	if (debug)
		std::cout << "Text: " << text << "(" << revText << "), Pattern: " << pattern << std::endl
				<< "FM1: " << fm1.vDesc.range.i1 << "-" << fm1.vDesc.range.i2 << ", FM2: " <<
					fm2.vDesc.range.i1 << "-" << fm2.vDesc.range.i2 << std::endl;

	bool res1 = true, res2 = true, res3 = true, res4 = true;

	for (i = 0; i < patternLength; ++i)
	{
		res1 &= goDown(fm1, pattern[patternLength - i - 1]);
		res2 &= goDown(fm2, pattern[i]);
		res3 &= rightExtend(bifm1, pattern[i]);
		res4 &= leftExtend(bifm2, pattern[patternLength - i - 1]);

		if (debug)
			std::cout << "Res: " << res1 << ", " << res2 << ", " << res3 << ", " << res4 << std::endl;

		if (debug)
			std::cout << "FM1: " << fm1.vDesc.range.i1 << "-" << fm1.vDesc.range.i2 << " (" << pattern[patternLength - i - 1] << "), FM2: " <<
			fm2.vDesc.range.i1 << "-" << fm2.vDesc.range.i2 << " (" << pattern[i] << ")" << std::endl << "- BFM1 fwd: " <<
			bifm1.fwdIter.vDesc.range.i1 << "-" << bifm1.fwdIter.vDesc.range.i2 << " (" << pattern[i] << "), BFM1 bwd: " <<
			bifm1.bwdIter.vDesc.range.i1 << "-" << bifm1.bwdIter.vDesc.range.i2 << " (" << pattern[i] << ")" << std::endl << "- BFM2 fwd: " <<
			bifm2.fwdIter.vDesc.range.i1 << "-" << bifm2.fwdIter.vDesc.range.i2 << " (" << pattern[patternLength - i - 1] << "), BFM2 bwd: " <<
			bifm2.bwdIter.vDesc.range.i1 << "-" << bifm2.bwdIter.vDesc.range.i2 << " (" << pattern[patternLength - i - 1] << ")" << std::endl;

		if (i == patternLength-1)
		{
			if (!(res1 == res2 && res2 == res3 && res3 == res4))
			{
				std::cerr << "Error1 (" << text << ", " << pattern << ")!" << std::endl;
				return 1;
			}
			else if (res1 == 1 && !( // all res are not 0!
				fm1.vDesc.range.i1 == bifm1.fwdIter.vDesc.range.i1 &&
				fm1.vDesc.range.i1 == bifm2.fwdIter.vDesc.range.i1 &&
				fm2.vDesc.range.i1 == bifm1.bwdIter.vDesc.range.i1 &&
				fm2.vDesc.range.i1 == bifm2.bwdIter.vDesc.range.i1 &&
				fm1.vDesc.range.i2 == bifm1.fwdIter.vDesc.range.i2 &&
				fm1.vDesc.range.i2 == bifm2.fwdIter.vDesc.range.i2 &&
				fm2.vDesc.range.i2 == bifm1.bwdIter.vDesc.range.i2 &&
				fm2.vDesc.range.i2 == bifm2.bwdIter.vDesc.range.i2
			))
			{
				std::cerr << "Error2 (" << text << ", " << pattern << ")!" << std::endl;
				return 1;
			}
		}
	}

	return 0;
}
template <typename TText, class TOccSpec, class TIndexSpec, typename TSpec>
inline bool
_rightExtendWithLargestRange(Iter<Index<TText, BidirectionalFMIndex<TOccSpec, TIndexSpec> >, VSTree<TopDown<ParentLinks<TSpec> > > > & it) {
	typedef typename Value<TText>::Type						TChar;
	typedef Index<TText, FMIndex<TOccSpec, TIndexSpec> >	TIndex;
	typedef typename Value<TIndex>::Type                	TAlphabet;

	TChar charWithLargestRange = 0;
	int largestRangeSize = 0;

	int alpSize = ValueSize<TAlphabet>::VALUE;
	for (int i = 0; i < alpSize; ++i)
	{
		TChar c = TAlphabet(i);

		int i1 = value(it.fwdIter).range.i1;
		int i2 = value(it.fwdIter).range.i2;

		if (!rightExtend(it, c))
			continue;

		if (i2-i1 > largestRangeSize)
		{
			charWithLargestRange = c;
			largestRangeSize = i2-i1;
		}

		goUp(it.fwdIter);
	}

	if (largestRangeSize == 0)
	{
		std::cout << "Feeeehler, wieso wurde denn überhaupt _rightExtendWithLargestRange aufgerufen?" << std::endl;
		return false;
	}

	std::cout << "extende mit " << charWithLargestRange << " (" << largestRangeSize << ")" << std::endl;
	rightExtend(it, charWithLargestRange);
	return true;
}

inline
double GetMedian(unsigned int daArray[], unsigned int iSize) {
    // Allocate an array of the same size and sort it.
	unsigned int* dpSorted = new unsigned int[iSize];
    for (unsigned int i = 0; i < iSize; ++i) {
        dpSorted[i] = daArray[i];
    }
    for (unsigned int i = iSize - 1; i > 0; --i) {
        for (unsigned int j = 0; j < i; ++j) {
            if (dpSorted[j] > dpSorted[j+1]) {
            	unsigned int dTemp = dpSorted[j];
                dpSorted[j] = dpSorted[j+1];
                dpSorted[j+1] = dTemp;
            }
        }
    }

    // Middle or average of middle values in the sorted array.
    unsigned int dMedian = 0.0;
    if ((iSize % 2) == 0) {
        dMedian = (dpSorted[iSize/2] + dpSorted[(iSize/2) - 1])/2.0;
    } else {
        dMedian = dpSorted[iSize/2];
    }
    delete [] dpSorted;
    return dMedian;
}

inline
double GetMode(unsigned int daArray[], unsigned int iSize) {
    // Allocate an int array of the same size to hold the
    // repetition count
	unsigned int* ipRepetition = new unsigned int[iSize];
    for (unsigned int i = 0; i < iSize; ++i) {
        ipRepetition[i] = 0;
        unsigned int j = 0;
        while ((j < i) && (daArray[i] != daArray[j])) {
            if (daArray[i] != daArray[j]) {
                ++j;
            }
        }
        ++(ipRepetition[j]);
    }
    unsigned int iMaxRepeat = 0;
    for (unsigned int i = 1; i < iSize; ++i) {
        if (ipRepetition[i] > ipRepetition[iMaxRepeat]) {
            iMaxRepeat = i;
        }
    }
    delete [] ipRepetition;
    return daArray[iMaxRepeat];
}

/*inline void
analyseUniformnessHumanGenome()
{
	typedef String<char> TText;
	typedef Index<TText, FMIndex<> > TIndex;
	typedef Iterator<TIndex, TopDown<ParentLinks<> > >::Type TIter;
	CharString id;
	TText read;
	TText genome;

	//const char *indexFile = "/home/chris/PhD/humanGenomeChrom1/indexChrom";

	std::cout << "[" <<  currentDateTime() << "]: Program startet." << '\n';

	// Read Chromosome
	SeqFileIn genomeFile("/media/chris/data/phd_data/human_g1k_v37.fasta");
	readRecord(id, genome, genomeFile);
	std::cout << "[" <<  currentDateTime() << "]: Chromosome loaded. Length = " << length(genome) << '\n';
	close(genomeFile);
	TIndex index2(genome);
	TIter iter2(index2);
	//indexRequire(index2, FibreSA());
	std::cout << "[" <<  currentDateTime() << "]: Index created." << '\n';
	//save(index, indexFile, OPEN_WRONLY | OPEN_CREATE);
	//std::cout << "[" <<  currentDateTime() << "]: Index saved to disk." << '\n';

	// Index Chromosome
	//genome = "";
	//TIndex index2(genome);
	////open(index2, indexFile);
	////std::cout << "[" <<  currentDateTime() << "]: Index loaded from disk." << '\n';
	//TIter iter2(index2);

	unsigned int const seedLength = 13;
	unsigned int const seedOffset = 1;

	//500 - seedLEngth + Offset = 491

	unsigned int _len = 0;

	SeqFileIn seqFileIn("/media/chris/data/phd_data/ERR843997.fastq");
	while (!atEnd(seqFileIn))
	{
		unsigned int seedsFound = 0;
		unsigned int seedsNotFound = 0;
		unsigned int numberHits[500-seedLength+seedOffset];
		readRecord(id, read, seqFileIn);
		// for each read
		for (unsigned int i = 0; i < length(read) - seedLength + 1; i += seedOffset)
		{
			// search seed (exact match)
			for (_len = 0; _len < seedLength; ++_len)
			{
				if(!goDown(iter2, read[i+_len])) {
					break;
				}
			}
			if (_len == seedLength)
			{
				numberHits[seedsFound] = countOccurrences(iter2);
				++seedsFound;
			}
			else
				++seedsNotFound;
			goRoot(iter2);
		}

		unsigned int minHits = 99999999;
		unsigned int maxHits = 0;
		unsigned int avgHits = 0;
		for (unsigned int i = 0; i < seedsFound; ++i)
		{
			if (numberHits[i] < minHits)
				minHits = numberHits[i];
			if (numberHits[i] > maxHits)
				maxHits = numberHits[i];
			avgHits += numberHits[i];
		}
		avgHits /= seedsFound;

		std::cout << "Found: " << seedsFound << "\tNot Found: " << seedsNotFound << "\tMin: " << minHits << "\tMax: " << maxHits << "\tAvg: " << avgHits <<
			"\tMedian: " << GetMedian(numberHits, seedsFound) << "\tMode: " << GetMode(numberHits, seedsFound) << std::endl;
	}
	close(seqFileIn);
	std::cout << "[" <<  currentDateTime() << "]: Program terminated." << '\n';
}*/

/*inline void
testRamUsage(bool const bidirectional = false)
{
	typedef String<char> TText;
	typedef Index<TText, FMIndex<> > TIndex;
	typedef Iterator<TIndex, TopDown<ParentLinks<> > >::Type TIter;
	CharString id;
	TText seq;

	const char *fileName = "/home/chris/indexTest/myindex";

	//std::cout << "[" <<  currentDateTime() << "]: Start" << std::endl;
	SeqFileIn seqFileIn("/home/chris/PhD/human_g1k_v37.fasta");
	readRecord(id, seq, seqFileIn);
	//std::cout << "[" <<  currentDateTime() << "]: Length of chromosome is " << length(seq) << '\n';
	close(seqFileIn);
	TIndex bifmIndex(seq);
	//std::cout << "[" <<  currentDateTime() << "]: Index done." << std::endl;

	TIter it(bifmIndex);
	//std::cout << "[" <<  currentDateTime() << "]: Iterator done." << std::endl;

	goDown(it, 'A');

	//delete seqFileIn;
	//delete *seq;

	//indexRequire(bifmIndex, FibreSA());

	//system("exec rm -r /home/chris/indexTest/ENTFERNE_MICH_WEGEN_KOMMENTAR_SYMBOL*");

	//save(bifmIndex, fileName);


	std::cout << "Index created" << std::endl;

	for (unsigned int i = 0; i < 100000; ++i) {
		goDown(it, 'A');
		//std::cout << "[" <<  currentDateTime() << "]:" << "Range:\t fwd: " << range(it.fwdIter).i1 << "-" << range(it.fwdIter).i2 << "\t bwd: " << range(it.bwdIter).i1 << "-" << range(it.bwdIter).i2 << std::endl;

		goDown(it, 'A');
		//std::cout << "[" <<  currentDateTime() << "]:" << "Range:\t fwd: " << range(it.fwdIter).i1 << "-" << range(it.fwdIter).i2 << "\t bwd: " << range(it.bwdIter).i1 << "-" << range(it.bwdIter).i2 << std::endl;
	}

	while(true) {}
}*/

template <typename TText, class TOccSpec, class TIndexSpec, typename TSpec>
inline void
find(
	Iter<Index<TText, BidirectionalFMIndex<TOccSpec, TIndexSpec> >, VSTree<TopDown<ParentLinks<TSpec> > > > &db,
	StringSet<TText> &sequences,
	unsigned int seedLength,
	unsigned int seedOffset,
	bool (*onAbortLambda)(unsigned int noErrors, unsigned int totalHits, unsigned int readLength, unsigned int foundReadLength, unsigned int currentHits),
	void (*onMatchLambda)(unsigned int pos, unsigned int noErrors)
) {
	typedef typename Iterator<StringSet<TText> >::Type TIterator;
	// for each sequence
	for (TIterator seqIter = begin(sequences, Standard()); seqIter != end(sequences, Standard()); ++seqIter)
	{
		unsigned int totalHits = 0;
		std::cout << "Sequence: " << *seqIter << std::endl;
		// for each seed
		for (unsigned int i = 0; i < length(*seqIter) - seedLength + 1; i += seedOffset)
		{
			unsigned int len = 0;
			unsigned int noErrors = 0;
			unsigned int startPos = i; // needed for going left

			std::cout << "\tSeed: " << infix(*seqIter, startPos, startPos + seedLength) << std::endl;

			// current seed has location [startPos, startPos+len-1]

			// search seed (exact match)
			for (len = 0; len < seedLength; ++len) {
				if(!rightExtend(db, (*seqIter)[i+len])) {
					std::cout << "\t\tSeed not found!" << std::endl;
					break;
				}
			}

			if (len != seedLength)
				continue; // skip this seed, since we did not find an exact match of it

			std::cout << "\t\tSeed found! Pos in Text: ";
			for (unsigned int j = 0; j < countOccurrences(db.fwdIter); ++j)
				std::cout << getOccurrences(db.fwdIter)[j] << ", ";
			std::cout << std::endl;

			unsigned int currentHits = countOccurrences(db.fwdIter);

			while (!(*onAbortLambda)(noErrors, totalHits, len, length(*seqIter), currentHits)) {
				currentHits = countOccurrences(db.fwdIter);

				//infix(str, begin, end): [begin..end-1], d.h. infix(str, 0, 3) gibt die ersten 3 positionen zurück: 0,1,2
				// direkter char-zugriff mit [i]: von 0 bis len-1

				// try to go right
				if (startPos+len < length(*seqIter) && rightExtend(db, (*seqIter)[startPos+len])) {
					++len;
					std::cout << "\t\tExtension: right (exact): " << infix(*seqIter, startPos, startPos + len);
				}
				// try to go left
				else if (startPos > 0 && leftExtend(db, (*seqIter)[startPos-1])) {
					--startPos;
					++len;
					std::cout << "\t\tExtension: left (exact): " << infix(*seqIter, startPos, startPos + len);
				}
				// approximate to the right
				else if (startPos+len < length(*seqIter)) {

					// So far: always take the first range
					//std::list<Pair<int> > ranges;
					// TODO: welchen character nehmen? alle durchprobieren? approx() wiederverwenden!
					//TText charStr = "";
					//appendValue(charStr, (*seqIter)[startPos+len]);
					//rightApproxExtend(db, 1, charStr, ranges); // TODO: +-1?
					// TODO: rechts/links angekommen -> abbrechen! da onAbortLambda nicht abbricht, werden die Treffer auch nicht gemeldet.
					_rightExtendWithLargestRange(db);
					++noErrors;
					++len;
					std::cout << "\t\tExtension: right (approx) " << infix(*seqIter, startPos, startPos + len);
				}
				else {
					currentHits = 0;
					std::cout << "\t\tExtension: not possible (" << startPos+len << " < " << length(*seqIter) - 1 << ")" << std::endl;
					break;
				}
				currentHits = countOccurrences(db.fwdIter);

				std::cout << " Pos in Text: ";
				for (unsigned int j = 0; j < currentHits; ++j)
					std::cout << getOccurrences(db.fwdIter)[j] << ", ";
				std::cout << std::endl;
			}

			std::cout << "\t\tExtension aborted with currentHits = " << currentHits << std::endl;

			if (currentHits > 0) {
				totalHits += currentHits;
				for (unsigned int j = 0; j < currentHits; ++j)
					onMatchLambda(getOccurrences(db.fwdIter)[j], noErrors);
			}

			goRoot(db);
			/*
			while (i+len < length(*seqIter) && rightExtend(db, (*seqIter)[i+len])) {
				++len;
				if (len >= seedLength) {
					unsigned int currentHits = countOccurrences(db.fwdIter);

					std::cout << infix(*seqIter, i, i + len) << " (" << range(db.fwdIter).i1 << "-" << range(db.fwdIter).i2 << "): ";
					for (unsigned j = 0; j < currentHits; ++j)
						std::cout << getOccurrences(db.fwdIter)[j] << ", ";
					std::cout << std::endl;

					if ((*onAbortLambda)(noErrors, totalHits, len, currentHits))
					{
						totalHits += currentHits;
						for (unsigned j = 0; j < currentHits; ++j)
							onMatchLambda(getOccurrences(db.fwdIter)[j], noErrors);
					}
				}
			}*/
		}
	}
}

}

#endif /* FIND_INDEX_TODO_H_ */
