// Method 1: brand new movieData. We the user clicks on 'finish':
// 1 ask the user where he wants to save the movie data (should provide a file
//   name and a path. Default filename = movieData.mat)
// 2 get the channelPaths, the pixelSize and the time interval from the GUI
// 3 create a movieData object (call the construction inside a try/catch
//   block statement since it can throw exception).
// 4 save it to disk: save([path filesep filename], 'movieData');
// 5 call the package controller with the current movieData object

// Method 2: from an existing movieData:
// 1 [path filename] = ask the user to select a movieData file
// 2 load([path filesep filename);
// 3 movieData.sanityCheck(path, filename, false); (again, call this method in a
//   try/catch block statement)
// 4 if no error, set the GUI components with the info found in movieData (GUI
//   components should be set to disabled).
// 5 when the user clicks on 'finish' call the package controller with the
//   current movieData object.

public class MovieData {
	
	public MovieData(String movieDataPath,
					 String movieDataFileName,
					 String[] channelPaths,
					 double pixelSize,
					 double timeInterval) throws Exception
	{
		movieDataPath_ = movieDataPath;
		movieDataFileName_ = movieDataFileName;
		nFrames_ = -1;
		imSize_[0] = -1;
		channelPaths_ = channelPaths;
		pixelSize_ = pixelSize;
		timeInterval_ = timeInterval;
		
		// To get the number of frames and the size of each image, we need first
		// to check that:
		// - each channel directory has the same number of frames
		// - all images have the same size (read only the tiff header for that)
		// These tests are performed by sanityCheck. the sanityCheck function
		// takes care of assigning imSize_ and nFrames_ if the tests are good.
		
		sanityCheck(movieDataPath, movieDataFileName);
	}

	// Get the path where the movieData is saved
	public String getMovieDataPath() {
		return movieDataPath_;
	}
	
	// Get the filename of the movieData
	public String getMovieDataFileName() {
		return movieDataFileName_;
	}
	
	// Get the number of images
	public int[] getImSize() {
		return imSize_;
	}
	
	// Get the number of frames
	public int getNFrames() {
		return nFrames_;
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
	// Any error encountered thrown an exception. It should be then used in a
	// try-catch statement. It begins by checking the validity of all data in
	// the movieData (movieDataPath_, imSize_, ...) and then, it calls the
	// sanityCheck function of each package in the package array. The 'full'
	// argument gives the ability to have 2 kind of check, one light and one
	// thorough.
	public void sanityCheck(String movieDataPath,
							String movieDataFileName,
							boolean full = false) throws Exception {
		// Check if the path and filename stored in the movieData are the same
		// than the ones provided in argument. They can differ if the movieData
		// file has been rename, move or copy to another location.
		if (movieDataPath_ != movieDataPath)
			movieDataPath_ = movieDataPath_;
		
		if (movieDataFileName_ != movieDataFilename)
			movieDataFileName_ = movieDataFileName;
		
		// Check that every channelPath exists
		// TO DO
		
		// Check that each channelPath contains the same number of image files.
		// During the check, get the actual number of frames on disk (nFrames).
		// TO DO

		if (nFrames_ != -1) {
			if (nFrames_ != nFrames)
				; // through exception: the number of frames stored in this
			      // movieData is different from the number of images in
			      // channel directories.
		} else {
			nFrames_ = nFrames;
		}
		
		// Check that each image in every channels has the same size. During the
		// check, get the size of the image.
		// TODO
		
		if (imSize_[0] != -1) {
			if (imSize_[0] != imSize[0] || imSize_[1] != imSize[1])
				; // through exception: the image size stored in this movieData
			      // is different from the image size of the image files.
		} else {
			imSize_ = imSize;
		}
		
		// Finally, call the sanity check on each package
		for (int i = 0; i < packages_.length(); ++i)
			packages_[i].sanityCheck(full);		
	}
	
	// This method may be called when a help button is triggered.
	public String getHelp() {
		return String("To set up a movieData, enter the pixel size and ...");
	}
	
	// Below starts the private field section
	
	private String movieDataPath_;
	private String movieDataFileName_;
	private int nFrames_;
	private int[] imSize_;
	private String[] channelPaths_;
	private double pixelSize_;
	private double timeInterval_;	
	private Process[] allProcesses_;
	private Package[] allPackages_;
}