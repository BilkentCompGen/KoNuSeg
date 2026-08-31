#include <mergeSegments.hpp>

// Function to insert a value into the array
void insert(std::vector<float>& arr, float value) {
    arr.push_back(value);
}

// Function to clear the array
void clear(std::vector<float>& arr) {
    arr.clear();
}




// Function to sort the array, eliminate the last 25% of the highest numbers,
// and return the median of the remaining values
float processData(std::vector<float>& arr) {
    // Sort the array
    std::sort(arr.begin(), arr.end());

    // Eliminate the last 25% of the highest numbers
    size_t size = arr.size();
    size_t remove_count = size / 4;  // Eliminate the top 25%
    size_t new_size = size - remove_count;

    // Create a new vector with the remaining values
    std::vector<float> remaining(arr.begin(), arr.begin() + new_size);

    // Find the median of the remaining values
    size_t remaining_size = remaining.size();
    if (remaining_size == 0) {
        std::cerr << "Error: No values to calculate median" << std::endl;
        return 0.0f; // Return a default value if no numbers are available
    }

    if (remaining_size % 2 == 0) {
        // If even, median is the average of the two middle values
        return roundCount((remaining[remaining_size / 2 - 1] + remaining[remaining_size / 2]) / 2.0f);
    } else {
        // If odd, median is the middle value
        return roundCount(remaining[remaining_size / 2]);
    }
}



void mergeSegments(
	string outfile_path, uint k_len, uint min_seg_len, uint count_threshold)
{
	ofstream wout(outfile_path + "_merged.bed");
	string currentSegment = "";
	uint segStartBase;
	uint segEndBase;
	float segCpNo;
	uint segCount;

	std::vector<float> values;
	uint coverage_depth = roundCount(static_cast<float>(k_reads) / k_assembly);
	uint gap_len = k_len;

		
	SEGMENTINFO *valSegment= new SEGMENTINFO();
	
	for (int pos = 0; pos < segTable->size(); pos++)
	{
		SEGMENTINFO* segItem = segTable->at(pos);
		// Each row in the cTable contains all the segments of the same contig

		for (int i = 0; i < segItem->size(); i++){
						
			// Initialising the first variables at the start of iteration through small segments
			if(currentSegment == "" )
			{
				currentSegment = segItem->at(i).segmentName;
				segStartBase = segItem->at(i).start;
				segEndBase = segItem->at(i).end;
				segCpNo = roundCount(segItem->at(i).count);	
				insert(values, segCpNo);

			}
			// If we are in the same contig and we have a non empty segment list
			else
			{ 
				if (currentSegment.compare(segItem->at(i).segmentName)==0)
				{	
					// If the segment at hand is smaller than the minimum requirement for the min size && the gap between the end of 
					// current segment and new one to be added is not > 3*k_len
					if ((segEndBase - segStartBase < min_seg_len) && (segItem->at(i).start - segEndBase)<= 3*k_len)
					{	
						segEndBase = segItem->at(i).end;
						insert(values, roundCount(segItem->at(i).count));
					}
					// If the segment is already >= than the minimum size requirement(MIN_SEG_SIZE)
					else
					{
						
						valSegment->push_back(SegmentInfo(currentSegment, segStartBase, segEndBase , processData(values)));
				
						currentSegment = segItem->at(i).segmentName;
						segStartBase = segItem->at(i).start;
						segEndBase = segItem->at(i).end;
						clear(values);
						insert(values, roundCount(segItem->at(i).count));		
					}
				}
				// If we are in a different contig
				else
				{		
					if (roundCount(segCpNo)>0)
					{	
						valSegment->push_back(SegmentInfo(currentSegment, segStartBase, segEndBase , processData(values)));
					}
					currentSegment = segItem->at(i).segmentName;
					segStartBase = segItem->at(i).start;
					segEndBase = segItem->at(i).end;
					clear(values);
					insert(values, roundCount(segItem->at(i).count));
				}
			}													
		}
	}

	//Inserting the last segment to the list
	valSegment->push_back(SegmentInfo(currentSegment, segStartBase, segEndBase , roundCount(segCpNo)));



	for (int j = 0; j < valSegment->size(); j++)
	{
		wout << valSegment->at(j).segmentName << "\t" <<valSegment->at(j).start<< "\t" <<valSegment->at(j).end<< "\t" <<valSegment->at(j).count<< endl;
	}

}

