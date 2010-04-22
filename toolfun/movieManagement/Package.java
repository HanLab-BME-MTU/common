// The first thing to do when a package control panel is opened is to check
// whether an object of the same class exists in the list of packages.
//
// If the package already exists:
// - a pop up windows should notify the user saying that the current package has
// already been used on this data set and the previous parameter set and results
// will be used.
// - the sanity check method of the package should be called. Any exception
// thrown during that check should set the flag icons (good,warning,wrong) next
// to each process.
//
// If the package does not exists:
// - call the constructor of the pacakge. Within the constructor of the class
// Package, we check whether a required process has not been already computed by
// another Package. If so, a pop window should ask the user whether he wants to
// use the already computed process or he prefer to create a new one. Then, we
// run the sanity c

// This class defines the abstract class Package from which every user-defined
// package will inherit.
public abstract class Package {
	// Constructor
	protected Package(movieData owner, String name)
	{
		owner_ = owner;
		name_ = name;
	}
	
	// Get the name of the package
	public String getName() {
		return name_;
	}
	
	// make a sanity check of the process
	public abstract void sanityCheck(boolean full);
	
	// more abstract methods...
	
	protected abstract Process[] getDefaultProcessList_();
	protected abstract Matrix getProcessDependenyMatrix_();
	
	// protected methods section
	protected void checkProcessDependencies_(boolean full) throws Exception {
	}
	
	// protected fields section
	protected movieData owner_;
	
	// private fields
	private String name_;
}