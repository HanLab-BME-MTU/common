// Here is an example of a concrete process (i.e. that implements the Process
// abstract class)
public class SegmentationProcess extends Process
{
	// Define the constructor
	public SegmentationProcess(movieData owner, String dateTime) {
		super(owner, "SegmentationProcess", dateTime)
	}
	
	// Get the ith mask path
	public String getMaskPath(int i) {
		return maskPaths_[i];
	}
	
	// Get the name of the function used to create the mask
	public String getFunctionName() {
		return functionName_;
	}
	
	// Get the list of parameters used in the function to create the mask
	public String getFunctionParams() {
		return functionParams_;
	}
	
	public void sanityCheck(boolean full) {
		// check that the maskPaths_ array is the same size that channelPaths_
		// array.
		
		if (maskPaths_.length() != owner_.channelPaths().length()) {
			// error
			
			// check mask path for every channel
			// check mask number for every channel == owner_.nFrames()
			
		}
	}
	
	// Private field section
	
	private String[] maskPaths_;
	private String functionName_;
	private String[] functionParams_;
	
}