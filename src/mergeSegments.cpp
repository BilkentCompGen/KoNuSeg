#include <mergeSegments.hpp>

// Function to insert a value into the array
void insert(std::vector<float>& arr, float value) {
    arr.push_back(value);
}

// Function to clear the array
void clear(std::vector<float>& arr) {
    arr.clear();
}

// Sort pooled copy estimates, discard the top quartile (at least 1 highest outlier when size >= 2)
// to suppress Bloom filter hash-collision inflation, and return the rounded median of the remainder.
float processData(std::vector<float>& arr) {
    if (arr.empty()) {
        std::cerr << "Error: No values to calculate median" << std::endl;
        return 0.0f;
    }

    std::sort(arr.begin(), arr.end());

    const size_t size = arr.size();
    const size_t remove_count = (size >= 2) ? ((size + 3) / 4) : 0;
    const size_t remaining_size = size - remove_count;

    if (remaining_size % 2 == 0) {
        return roundCount((arr[remaining_size / 2 - 1] + arr[remaining_size / 2]) / 2.0f);
    } else {
        return roundCount(arr[remaining_size / 2]);
    }
}


bool mergeSegments(
	string outfile_path, uint k_len, uint min_seg_len, uint count_threshold)
{
	string merged_path = outfile_path + "_merged.bed";
	ofstream wout(merged_path);
	if (!wout.is_open()) {
		cerr << "[!] Error: Could not open merged output BED file: " << merged_path << endl;
		return false;
	}

	string currentSegment = "";
	uint segStartBase = 0;
	uint segEndBase = 0;

	std::vector<float> values;
	SEGMENTINFO valSegment;
	
	for (size_t pos = 0; pos < segTable->size(); pos++)
	{
		SEGMENTINFO* segItem = segTable->at(pos);
		if (segItem == nullptr) continue;

		// Each row in segTable contains all the segments of the same contig
		for (size_t i = 0; i < segItem->size(); i++)
		{
			const SegmentInfo& nextSeg = segItem->at(i);
			const float nextRounded = roundCount(nextSeg.count);

			if (currentSegment.empty())
			{
				currentSegment = nextSeg.segmentName;
				segStartBase = nextSeg.start;
				segEndBase = nextSeg.end;
				insert(values, nextRounded);
			}
			else if (currentSegment == nextSeg.segmentName)
			{
				const uint currentLen = (segEndBase >= segStartBase) ? (segEndBase - segStartBase) : 0;
				const uint nextLen = (nextSeg.end >= nextSeg.start) ? (nextSeg.end - nextSeg.start) : 0;
				const bool withinGap = (nextSeg.start <= segEndBase) ||
					((nextSeg.start - segEndBase) <= 3 * k_len);

				// Avoid absorbing a >= k_len flanking locus into an already full-sized (>= min_seg_len) block
				// when their copy numbers differ by more than count_threshold, while smoothing sub-kmer (< k_len) artifacts.
				const bool incompatibleLongBlock = (currentLen >= k_len) &&
					(nextLen >= min_seg_len) &&
					!values.empty() &&
					(std::fabs(nextRounded - values.back()) > static_cast<float>(count_threshold));

				if (currentLen < min_seg_len && withinGap && !incompatibleLongBlock)
				{
					if (nextSeg.end > segEndBase) {
						segEndBase = nextSeg.end;
					}
					insert(values, nextRounded);
				}
				else
				{
					float mergedCount = processData(values);
					if (mergedCount > 0.0f) {
						valSegment.push_back(SegmentInfo(currentSegment, segStartBase, segEndBase, mergedCount));
					}

					currentSegment = nextSeg.segmentName;
					segStartBase = nextSeg.start;
					segEndBase = nextSeg.end;
					clear(values);
					insert(values, nextRounded);
				}
			}
			// Different contig
			else
			{
				float mergedCount = processData(values);
				if (mergedCount > 0.0f)
				{
					valSegment.push_back(SegmentInfo(currentSegment, segStartBase, segEndBase, mergedCount));
				}
				currentSegment = nextSeg.segmentName;
				segStartBase = nextSeg.start;
				segEndBase = nextSeg.end;
				clear(values);
				insert(values, nextRounded);
			}
		}
	}

	// Flush the last segment using processData(values) when non-empty
	if (!currentSegment.empty() && !values.empty())
	{
		float mergedCount = processData(values);
		if (mergedCount > 0.0f) {
			valSegment.push_back(SegmentInfo(currentSegment, segStartBase, segEndBase, mergedCount));
		}
	}

	for (size_t j = 0; j < valSegment.size(); j++)
	{
		wout << valSegment.at(j).segmentName << "\t" << valSegment.at(j).start << "\t" << valSegment.at(j).end << "\t" << valSegment.at(j).count << "\n";
	}

	wout.flush();
	bool ok = wout.good();
	wout.close();
	return ok;
}
