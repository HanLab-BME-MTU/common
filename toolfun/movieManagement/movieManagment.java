//
//  movieManagment.java
//  
//
//  Created by Sylvain Berlemont on 3/22/10.
//  Copyright 2010 Danuser Lab - Harvard Medical School. All rights reserved.
//

// This defines the abstract class Process from which every user-defined process
// will inherit.
public abstract class Process {
	// Contructor
	protected Process(movieManagement owner, String name)
	{
		owner_ = owner;
	}
		
	// Get the name of the process
	public String getName() {
		return name_;
	}
	
	// make a sanity check of the process
	public abstract boolean sanityCheck();
	
	// more abstract methods...
	
	// protected field section
	protected const movieManagement owner_;
	
	// private fields
	private String name_;
}

// Here is an example of a concrete process (i.e. that implements the Process abstract class)
public class maskProcess extends Process
{
	// Define the constructor
	public maskProcess(movieManagement owner) {
		super(owner, 'maskProcess')
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
		
		if maskPaths_.length() != owner_.channelPaths()
			// error
		
		// check mask path for every channel
	}
	
	// Private field section
	
	String[] maskPaths_;
	String functionName_;
	String[] functionParams_;
}

// This defines the abstract class Process from which every user-defined process
// will inherit.
public abstract class Package {
	// Contructor
	protected Package(movieManagement owner, String name)
	{
		owner_ = owner;
	}
	
	// Get the name of the package
	public String getName() {
		return name_;
	}
	
	// make a sanity check of the process
	public abstract boolean sanityCheck();
	
	// more abstract methods...
	
	// protected field section
	protected const movieManagement owner_;
	
	// private fields
	private String name_;
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
	
	// TODO: Define constructors (that might be specific to Matlab OO)

	// Return the number of frames
	public int nFrames() {
		return nFrames_;
	}

	// Get the current spatial unit (nm, um, ...)
	public String getSpatialUnit() {
		return spatialUnit_;
	}
	
	// Set the spatial unit.
	public void setSpatialUnit(String spatialUnit) {
		spatialUnit_ = spatialUnit;		
	}
	
	// Get the current temporal unit (s, min, h, ...)
	public String getTemporalUnit() {
		return temporalUnit_;
	}
	
	// Set the temporal unit.
	public void setTemporalUnit(String temporalUnit) {
		temporalUnit_ = temporalUnit;
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

	// Add a path to the channelPaths array
	public void addChannelPath(String channelPath) {
	}
	
	// Get the number of processes
	public int getNumberOfProcesses() {
		return processes_.length();
	}
	
	// Get the ith process
	public Process getProcess(int i) {
		return processes_[i];
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
	
	// Add a package to the package array
	public void addPackage(Packge p) {
	}
	
	// This method checks the validity of the data stored in that object.
	public boolean sanityCheck(boolean full) {
		boolean isValid = true;

		// check whether channelPaths exist and that each channelPath contains
		// nFrames tif files, ...
		
		// Call the sanity check of every process
		for (int i = 0; i < processes_.length(); ++i)
			isValid &= processes_[i].sanityCheck();
		
		if (~isValid)
			return false;
		
		// Call the sanity check of every package
		for (int i = 0; i < packages_.length(); ++i)
			isValid &= packages_[i].sanityCheck();
		
		return isValid;
	}
	
	// Below starts the private field section
	
	private int nFrames_;
	private String spatialUnit_;
	private String temporalUnit_;
	private double pixelSize_;
	private double timeInterval_;
	private String[] channelPaths_;
	
	private Process[] processes_;
	private Package[] pachages_;
}
