

// This defines the abstract class Process from which every user-defined process
// will inherit.
public abstract class Process {
	// Constructor
	protected Process(movieManagement owner, String name, String dateTime)
	{
		owner_ = owner;
	}
		
	// Get the name of the process
	public String getName() {
		return name_;
	}

	public String getDateTime() {
		return dateTime_;
	}		

	// make a sanity check of the process
	public abstract boolean sanityCheck();
	
	// more abstract methods...
	
	// protected field section
	protected const movieManagement owner_;
	protected String dateTime_;	

	// private fields
	private String name_;

}

// Here is an example of a concrete process (i.e. that implements the Process abstract class)
public class maskProcess extends Process
{
	// Define the constructor
	public maskProcess(movieManagement owner, String dateTime) {
		super(owner, 'maskProcess', dateTime)
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
	
	public boolean sanityCheck() {
		// check that the maskPaths_ array is the same size that channelPaths_
		// array.
		
		if maskPaths_.length() != owner_.channelPaths().length()
			// error
		
		// check mask path for every channel
		// check mask number for every channel == owner_.nFrames()

	}
	
	// Private field section
	
	String[] maskPaths_;
	String functionName_;
	String[] functionParams_;

}

// This defines the abstract class Process from which every user-defined process
// will inherit.
public abstract class Package {
	// Constructor
	protected Package(movieManagement owner, String name)
	{
		owner_ = owner;
	}
	
	// Get the name of the package
	public String getName() {
		return name_;
	}

	public String[] getProcesses() {
		return processes_;
	}
	
	// make a sanity check of the process
	public abstract boolean sanityCheck();
	
	// more abstract methods...
	
	// protected field section
	protected const movieManagement owner_;
	
	// private fields
	private String name_;
	private String[] processes_;

}

// Here is an example of a concrete package (i.e. that implements the Package abstract class)
public class bioSensorPackage extends Package {
	// Define the constructor
	public bioSensorPackage(movieManagement owner) {
		super(owner, 'bioSensorPackage')
	}
	
	public boolean sanityCheck() {
		// check that every process related to biosensor has been processed:
		
		// Traverse the list
		for (int i = 1; i < owner_.getNumberOfProcesses(); ++i) {
			if (owner_.getProcess(i).getName() != 'maskProcess') {				
				// error => Mask process hasn't been performed.
			}
			else {
				// maskProcess is here and needs to be sanity checked;
				// owner_.getProcess(i).sanityCheck()
			}
		}
	}
}
	

public class movieManagment {
	
	public movieManagement(int[] imSize,int nFrames,String[] channelPaths,double pixelSize,double timeInterval)

	// Return the number of frames
	public int nFrames() {
		return nFrames_;
	}

	//Return the number of images
	public int[] imSize() {
		return imSize_;
	}
	
	// Get pixel size
	public double getPixelSize() {
		return pixelSize_;
	}
	
	// Get time interval
	public double getTimeInterval() {
		return timeInterval_;
	}
	
	// Get the ith channel path
	public String getChannelPath(int i) {
		return channelPaths_[i];
	}
	
	// Get the number of processes
	public int getNumberOfProcesses() {
		return processes_.length();
	}
	
	// Get the ith process
	public Process getProcess(int i) {
		return processes_[i];
	}

	// Set the ith process
	public Process setProcess(int i, Process process) {
		processes_[i] = process;
	}
	
	// Add a process to the processes array
	public void addProcess(Process p) {
	}
	
	// Get the number of packages
	public int getNumberOfPackages() {
		reutrn packages_.length();
	}
	
	// Get the ith package
	public Package getPackage(int i) {
		return packages_[i];
	}
	
	//Set the ith package
	public Package setPackage(int i, Package package){
		packages_[i] = package;
	}
	
	// Add a package to the package array
	public void addPackage(Package p) {
	}
	
	// This method checks the validity of the data stored in that object.
	public boolean sanityCheck(boolean full) {
		boolean isValid = true;

		// check whether channelPaths exist and that each channelPath contains
		// nFrames tif files, ...		
		
		// Call the sanity check of every process
		for (int i = 0; i < processes_.length(); ++i)
			isValid &= processes_[i].sanityCheck(full);
		
		if (~isValid)
			return false;
				
		// Call the sanity check of every package
		for (int i = 0; i < packages_.length(); ++i)
			isValid &= packages_[i].sanityCheck(full);
		
		return isValid;
	}
	
	// Below starts the private field section
	
	private int nFrames_;
	private int[] imSize_;
	private String[] channelPaths_;
	private double pixelSize_;
	private double timeInterval_;	
	private Process[] processes_;
	private Package[] packages_;
}
