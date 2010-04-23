// The first thing to do when a package control panel is opened is to check
// whether an object of the same class exists in the list of packages.
//
// If the package already exists:
// - a pop up windows should notify the user that the current package has already
// been used on this data. The pop up windows should have 2 buttons 'Choose...'
// and 'Create'. 'Choose' will open another window with the list of existing related
// packages. 'Create' will create a new package.
// - the sanity check method of the package should be called. Any exception
// thrown during that check should set the flag icons (good,warning,wrong) next
// to each process.
//
// If the package does not exists:
// - call the constructor of the pacakge. Within the constructor of the class
// Package, we check whether a required process has not been already computed by
// another Package. If so, a pop window should ask the user whether he wants to
// use the already computed process or he prefers to create a new one. Again, this
// pop up window should have 2 buttons 'Choose...' and 'Create'. 'Choose' will
// open another window with the list of existing related processes. 'Create' will
// create a new process with default parameter.
// - A sanity check is then performed to check the consistency of the object. Any
// exception thrown during that check should set the flag icons (good,warning,wrong)
// next to each process.

// This class defines the abstract class Package from which every user-defined
// package will inherit.
public abstract class Package {
	// Constructor
	protected Package(movieData owner, String name)
	{
		owner_ = owner;
		name_ = name;
		
		// initialize localProcesses_ with the list of required process. Each
		// process are built with default parameters.
		localProcesses_ = getDefaultProcessList_();
		
		int nLocalProcesses = length(localProcesses_);
		int nGlobalProcesses = length(owner_.processes_);
		
		boolean[] isSame = new int[nProcesses];
		
		// Iterate over every required process and check whether processes of
		// the same class exists
		for (int i = 0; i < nLocalProcesses; ++i) {
			
			String className = className(localProcesses_[i]);
			
			for (int j = 0; j < nGlobalProcesses; ++j) {
				isSame = classname(owner_.processes_[j]) == className;
			}
			
			if (any(isSame)) {
				// At least one process is of the same class as the required ith
				// process. This is when we need to show the popup window with
				// "choose..." and "create" button (see header at the top of the
				// file). Every localProcesses_[j] where isSame[j] == 1 needs to
				// be displayed as a potential choice. 
				
				// TODO: ask user
				
				// test the user answer
				if (response == "choose..." && chosenProcess ~= end) {
					// Replace the default process by the one chosen by user
					localProcesses_[i] = owner_.processes_[chosenProcess];
				} else {
					// Add the new process into owner's process list. 
					owner_.processes_.add(localProcesses_[i]);
				}
			}
		}
		
		// Call the sanity check
		sanityCheck(false);
	}
	
	// Get the name of the package
	public String getName() {
		return name_;
	}
	
	// make a sanity check of the package
	public abstract void sanityCheck(boolean full);
	
	// get the list of processes required by the package. These processes contains
	// the default parameters to be used in the algorithm.
	protected abstract Process[] getDefaultProcessList_();
	
	// This method return the matrix of process dependency. It represents the
	// directed graph of processes. Only the upper triangular part of the matrix
	// is considered. If there is no dependency between processes, return a
	// zero matrix, where the size is equal to the number of processes. Any
	// process should have selft dependency, i.e. there should not be any
	// mat(i,i) == 1. The type of the matrix should be logical, i.e.
	// mat = false(nProcess); mat(i,j) = true, etc.
	protected abstract Matrix getProcessDependenyMatrix_();
	
	// This method traverse the list of local processes and call for each process
	// its sanityCheck methods. The list of processes is traversed using a  
	protected void checkProcessDependencies_(boolean full) throws Exception {
		
		Matrix depMatrix = getProcessDependenyMatrix_();
		
	}

	// private fields
	private String name_;
	
	// protected fields section
	protected movieData owner_;
	// list of required processes
	protected Process[] localProcesses_;
}